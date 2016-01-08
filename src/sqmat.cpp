#include "sqmat.hpp"
#include "lapack.hpp"

template<typename T>
SqMat<T>::SqMat(const Vector<T> &vec):
  Matrix<T>(vec.size(), vec.size())
{
  for(size_t idx = 0; idx < vec.size(); ++idx)
    this->at(idx, idx) = vec.at(idx);
}

template<typename T>
SqMat<T> &SqMat<T>::operator=(const SqMat<T> &rsqmat)
{
  Matrix<T>::nrow = rsqmat.size();
  Matrix<T>::ncol = rsqmat.size();
  Matrix<T>::vec  = rsqmat.Matrix<T>::vec; 
  
  return *this;
}

template<typename T>
SqMat<T> &SqMat<T>::operator=(const Matrix<T> &rmat)
{
  assert(rmat.row_size() == rmat.col_size());
  Matrix<T>::operator=(rmat);
  
  return *this;
}

template<typename T>
T SqMat<T>::trace() const
{
  T tr = T();
  
  for(size_t i = 0; i < this->size(); ++i)
    tr += this->at(i,i);
  
  return tr;
}

template<>
const SqMat<dreal> &SqMat<dreal>::inverse()
{
  lapack::dinvs(this->size(), this->data());
  
  return *this;
}

template<>
const SqMat<dcmplx> &SqMat<dcmplx>::inverse()
{
  lapack::zinvs(this->size(), this->data());
  
  return *this;
}

template<>
void SqMat<dreal>::diagonalize(SqMat<dreal> &sqmat, Vector<dreal> &eigval)
{
  eigval.resize(sqmat.size());
  
  lapack::dsyevd(sqmat.size(), sqmat.data(), eigval.data());
}

template<>
void SqMat<dreal>::diagonalize(size_t il, size_t ir, SqMat<dreal> &sqmat, Vector<dreal> &eigval)
{
  assert(il <= ir);
  assert(ir < sqmat.size());
  
  size_t m = ir - il + 1;
  
  eigval.resize(m);
  
  SqMat<dreal> eigvec(sqmat.size(), m);
  
  lapack::dsyevr(sqmat.size(), 0, 0, il+1, ir+1, sqmat.data(), eigval.data(), eigvec.data());
  
  for(size_t i = il; i <= ir; ++i)
    for(size_t j = 0; j < sqmat.size(); ++j)
      sqmat.at(j, i) = eigvec.at(j, i);
}

template<>
void SqMat<dcmplx>::diagonalize(SqMat<dcmplx> &sqmat, Vector<dreal> &eigval)
{
  eigval.resize(sqmat.size());
  
  lapack::zheevd(sqmat.size(), sqmat.data(), eigval.data());
}

template<>
void SqMat<dcmplx>::diagonalize(size_t il, size_t ir, SqMat<dcmplx> &sqmat, Vector<dreal> &eigval)
{
  assert(il <= ir);
  assert(ir < sqmat.size());
  
  size_t m = ir - il + 1;
  
  eigval.resize(m);
  
  SqMat<dcmplx> eigvec(sqmat.size(), m);
  
  lapack::zheevr(sqmat.size(), 0, 0, il+1, ir+1, sqmat.data(), eigval.data(), eigvec.data());
  
  for(size_t i = il; i <= ir; ++i)
    for(size_t j = 0; j < sqmat.size(); ++j)
      sqmat.at(j, i) = eigvec.at(j, i);
}

template class SqMat<dreal>;
template class SqMat<dcmplx>;
