#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include "precision.hpp"
#include "vector.hpp"

template<typename T>
class Matrix {
  
public:
  
  Matrix(): nrow(0), ncol(0), vec() {}
  
  Matrix(size_t nr, size_t nc, T value = T()):
    nrow(nr), ncol(nc), vec(nr*nc, value)
  {}
  
  virtual ~Matrix() {}
  
  size_t row_size() const
  { return nrow; }
  
  size_t col_size() const
  { return ncol; }
  
  T *data()
  { return vec.data(); }
  
  const T *data() const
  { return vec.data(); }
  
  T &at(size_t ir, size_t ic);
  
  const T &at(size_t ir, size_t ic) const;
  
  void resize(size_t nr, size_t nc, T value = T());
  
  Matrix<T> &operator=(const Matrix<T> &rmat);
  Matrix<T> operator+(const Matrix<T> &rmat) const;
  Matrix<T> operator-(const Matrix<T> &rmat) const;
  Matrix<T> operator*(const Matrix<T> &rmat) const;
  Vector<T> operator*(const Vector<T> &rvec) const;
  Vector<T> operator[](size_t ic) const;
  
  Matrix<T> transpose() const;
  
  static void vectorize(const Matrix<T> &mat, std::vector< Vector<T> > &vecstore);
  
  static void vectorize(size_t ibeg, size_t iend, const Matrix<T> &mat,
                        std::vector< Vector<T> > &vecstore);
  
protected:
  size_t    nrow;
  size_t    ncol;
  Vector<T> vec;
};

template<typename T>
std::ostream &operator<<(std::ostream &out, const Matrix<T> &mat);

typedef Matrix<dreal>  MatDreal;
typedef Matrix<dcmplx> MatDcmplx;

#endif
