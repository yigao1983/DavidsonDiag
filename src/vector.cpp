#include <iostream>
#include <cassert>
#include "vector.hpp"

template<typename T>
Vector<T> &Vector<T>::operator=(const Vector<T> &rvec)
{
  this->vector = rvec.vector;
  
  return *this;
}

template<typename T>
Vector<T> Vector<T>::operator+(const T &add) const
{
  Vector<T> vec = *this;
  
  for(size_t i = 0; i < vec.size(); ++i)
    vec.at(i) += add;
  
  return vec;
}

template<typename T>
Vector<T> Vector<T>::operator-(const T &sub) const
{
  Vector<T> vec = *this;
  
  for(size_t i = 0; i < vec.size(); ++i)
    vec.at(i) -= sub;
  
  return vec;
}

template<typename T>
Vector<T> Vector<T>::operator*(const T &fac) const
{
  Vector<T> vec = *this;
  
  for(size_t i = 0; i < vec.size(); ++i)
    vec.at(i) *= fac;
  
  return vec;
}

template<typename T>
Vector<T> Vector<T>::operator+(const Vector<T> &rvec) const
{
  assert(this->size() == rvec.size());
  
  Vector<T> vec = *this;
  
  for(size_t i = 0; i < vec.size(); ++i)
    vec.at(i) += rvec.at(i);
  
  return vec;
}

template<typename T>
Vector<T> Vector<T>::operator-(const Vector<T> &rvec) const
{
  assert(this->size() == rvec.size());
  
  Vector<T> vec = *this;
  
  for(size_t i = 0; i < vec.size(); ++i)
    vec.at(i) -= rvec.at(i);
  
  return vec;
}

template<typename T>
T Vector<T>::operator*(const Vector<T> &rvec) const
{
  assert(this->size() == rvec.size());
  
  T dp = T();
  
  for(size_t i = 0; i < this->size(); ++i)
    dp += this->at(i) * rvec.at(i);
  
  return dp;
}

template<>
dcmplx Vector<dcmplx>::operator*(const Vector<dcmplx> &rvec) const
{
  assert(this->size() == rvec.size());
  
  dcmplx dp = dcmplx();
  
  for(size_t i = 0; i < this->size(); ++i)
    dp += std::conj( this->at(i) ) * rvec.at(i);
  
  return dp;
}

template<typename T>
dreal Vector<T>::norm() const
{
  dcmplx sq = (*this) * (*this);
  
  return std::sqrt(sq.real());
}

template<typename T>
Vector<T> Vector<T>::conj() const
{
  return *this;
}

template<>
Vector<dcmplx> Vector<dcmplx>::conj() const
{
  Vector<dcmplx> vec = *this;
  
  for(size_t i = 0; i < vec.size(); ++i)
    vec.at(i) = std::conj(this->at(i)); 
  
  return vec;
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const Vector<T> &vec)
{
  for(size_t i = 0; i < vec.size(); ++i)
    out << vec.at(i) << " ";
  
  out << std::endl;
  
  return out;
}

template class Vector<dreal>;
template class Vector<dcmplx>;
template std::ostream &operator<<(std::ostream &out, const Vector<dreal> &vec);
template std::ostream &operator<<(std::ostream &out, const Vector<dcmplx> &vec);
