#ifndef SQMAT_HPP
#define SQMAT_HPP

#include <cassert>
#include "matrix.hpp"

template<typename T>
class SqMat : public Matrix<T> {
  
public:
  
  static void diagonalize(SqMat<T> &sqmat, Vector<dreal> &eigval);
  
  static void diagonalize(size_t il, size_t ir, SqMat<T> &sqmat, Vector<dreal> &eigval);
  
  SqMat(): Matrix<T>() {}
  
  SqMat(size_t nd, T value = T()):
    Matrix<T>(nd, nd, value)
  {}
  
  SqMat(const Matrix<T> &rmat)
  { 
    assert(rmat.row_size() == rmat.col_size());
    Matrix<T>::operator=(rmat);
  }
  
  SqMat(const Vector<T> &rvec);
  
  virtual ~SqMat() {}
  
  size_t size() const
  { return Matrix<T>::row_size(); }
  
  void resize(size_t nr, size_t nc, T value = T())
  {
    assert(nr == nc);
    resize(nr, value);
  }
  
  void resize(size_t nd, T value = T())
  { Matrix<T>::resize(nd, nd, value); }
  
  SqMat<T> &operator=(const SqMat<T> &rsqmat);
  SqMat<T> &operator=(const Matrix<T> &rmat);
  
  T trace() const;
  
  const SqMat<T> &inverse();
};

typedef SqMat<dreal>  SqMatDreal;
typedef SqMat<dcmplx> SqMatDcmplx;

#endif
