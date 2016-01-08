#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <vector>
#include "precision.hpp"

template<typename T>
class Vector {
  
public:
  Vector(): vector() {}
  
  Vector(size_t n, T value = T()): vector(n, value) {}
  
  virtual ~Vector() {}
  
  T &at(size_t idx)
  { return vector.at(idx); }
  
  const T &at(size_t idx) const
  { return vector.at(idx); }
  
  void clear()
  { vector.clear(); }
  
  T *data()
  { return vector.data(); }
  
  const T *data() const
  { return vector.data(); }
  
  void pop_back()
  { vector.pop_back(); }
  
  void push_back(const T &value)
  { vector.push_back(value); }
  
  void resize(size_t n, T value = T())
  { vector.resize(n, value); }
  
  size_t size() const
  { return vector.size(); }
  
  T &operator[](size_t idx)
  { return vector[idx]; }
  
  const T &operator[](size_t idx) const
  { return vector[idx]; }
  
  Vector<T> &operator=(const Vector<T> &rvec);
  Vector<T> operator+(const T &add) const;
  Vector<T> operator-(const T &sub) const;
  Vector<T> operator*(const T &fac) const;
  Vector<T> operator+(const Vector<T> &rvec) const;
  Vector<T> operator-(const Vector<T> &rvec) const;
  T operator*(const Vector<T> &rvec) const;
  dreal norm() const;
  Vector<T> conj() const;
  
protected:
  std::vector<T> vector;
};

template<typename T>
std::ostream &operator<<(std::ostream &out, const Vector<T> &vec);

typedef Vector<dreal>  VectorDreal;
typedef Vector<dcmplx> VectorDcmplx;

#endif
