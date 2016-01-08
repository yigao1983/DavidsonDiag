#include <cassert>
#include "matrix.hpp"
#include "lapack.hpp"

template<typename T>
T &Matrix<T>::at(size_t ir, size_t ic)
{
  assert(ir < nrow);
  assert(ic < ncol);
  
  return vec.at(ic*nrow + ir);
}

template<typename T>
const T &Matrix<T>::at(size_t ir, size_t ic) const
{
  assert(ir < nrow);
  assert(ic < ncol);
  
  return vec.at(ic*nrow + ir);
}

template<typename T>
void Matrix<T>::resize(size_t nr, size_t nc, T value)
{
  nrow = nr;
  ncol = nc;
  vec.resize(nr*nc, value);
}

template<typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &rmat)
{
  nrow = rmat.nrow;
  ncol = rmat.ncol;
  vec  = rmat.vec;
  
  return *this;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &rmat) const
{
  assert(nrow == rmat.nrow && ncol == rmat.ncol);
  
  Matrix<T> matsum = *this;
  
  for(size_t idx = 0; idx < matsum.vec.size(); ++idx)
    matsum.vec.at(idx) += rmat.vec.at(idx);
  
  return matsum;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &rmat) const
{
  assert(nrow == rmat.nrow && ncol == rmat.ncol);
  
  Matrix<T> matsub = *this;
  
  for(size_t i = 0; i < matsub.vec.size(); ++i)
    matsub.vec.at(i) -= rmat.vec.at(i);
  
  return matsub;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &rmat) const
{
  assert(ncol == rmat.nrow);
  
  Matrix<T> matprod(nrow, rmat.ncol);
  
  for(size_t ir = 0; ir < matprod.nrow; ++ir)
    for(size_t ic = 0; ic < matprod.ncol; ++ic)
      for(size_t j = 0; j < ncol; ++j)
        matprod.at(ir, ic) += this->at(ir, j) * rmat.at(j, ic);
  
  return matprod;
}

template<>
Matrix<dreal> Matrix<dreal>::operator*(const Matrix<dreal> &rmat) const
{
  assert(ncol == rmat.nrow);
  
  Matrix<dreal> matprod(nrow, rmat.ncol);
  
  lapack::dgemm(nrow, rmat.ncol, ncol, 'N', 'N', 1.0, this->data(), rmat.data(), 0.0, matprod.data());
  
  return matprod;
}

template<>
Matrix<dcmplx> Matrix<dcmplx>::operator*(const Matrix<dcmplx> &rmat) const
{
  assert(ncol == rmat.nrow);
  
  Matrix<dcmplx> matprod(nrow, rmat.ncol);
  
  lapack::zgemm(nrow, rmat.ncol, ncol, 'N', 'N', 1.0, this->data(), rmat.data(), 0.0, matprod.data());
  
  return matprod;
}

template<typename T>
Vector<T> Matrix<T>::operator*(const Vector<T> &rvec) const
{
  assert(ncol == rvec.size());
  
  Vector<T> vecprod(nrow);
  
  for(size_t ir = 0; ir < nrow; ++ir)
    for(size_t ic = 0; ic < ncol; ++ic)
      vecprod.at(ir) += this->at(ir, ic) * rvec.at(ic);
  
  return vecprod;
}

template<typename T>
Vector<T> Matrix<T>::operator[](size_t ic) const
{
  assert(ic < ncol);
  
  Vector<T> vecslice(nrow);
  
  for(size_t idx = 0; idx < nrow; ++idx)
    vecslice[idx] = this->at(idx, ic);
  
  return vecslice;
}

template<>
Vector<dreal> Matrix<dreal>::operator*(const Vector<dreal> &rvec) const
{
  assert(ncol == rvec.size());
  
  Vector<dreal> vecprod(nrow);
  
  lapack::dgemm(nrow, 1, ncol, 'N', 'N', 1.0, this->data(), rvec.data(), 0.0, vecprod.data());
  
  return vecprod;
}

template<>
Vector<dcmplx> Matrix<dcmplx>::operator*(const Vector<dcmplx> &rvec) const
{
  assert(ncol == rvec.size());
  
  Vector<dcmplx> vecprod(nrow);
  
  lapack::zgemm(nrow, 1, ncol, 'N', 'N', 1.0, this->data(), rvec.data(), 0.0, vecprod.data());
  
  return vecprod;
}

template<typename T>
Matrix<T> Matrix<T>::transpose() const
{
  Matrix<T> matt(ncol, nrow);
  
  for(size_t ir = 0; ir < nrow; ++ir)
    for(size_t ic = 0; ic < ncol; ++ic)
      matt.at(ic, ir) = this->at(ir, ic);
  
  return matt;
}

template<typename T>
void Matrix<T>::vectorize(const Matrix<T> &mat, std::vector< Vector<T> > &vecstore)
{
  vecstore.resize(mat.col_size());
  
  for(size_t ic = 0; ic < mat.col_size(); ++ic)
    vecstore[ic] = mat[ic];
}

template<typename T>
void Matrix<T>::vectorize(size_t ibeg, size_t iend, const Matrix<T> &mat,
                          std::vector< Vector<T> > &vecstore)
{
  assert(ibeg <= iend);
  assert(iend < mat.col_size());
  
  size_t m = iend - ibeg + 1;
  
  vecstore.resize(m);
  
  for(size_t i = 0; i < m; ++i) {
    size_t ic = ibeg + i;
    vecstore[i] = mat[ic];
  }
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const Matrix<T> &mat)
{
  for(size_t ir = 0; ir < mat.row_size(); ++ir) {
    for(size_t ic = 0; ic < mat.col_size(); ++ic)
      out << mat.at(ir, ic) << " ";
    out << std::endl;
  }
  
  return out;
}

template class Matrix<dreal>;
template class Matrix<dcmplx>;

template std::ostream &operator<<(std::ostream &out, const Matrix<dreal> &mat);
template std::ostream &operator<<(std::ostream &out, const Matrix<dcmplx> &mat);
