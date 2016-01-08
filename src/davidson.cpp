#include <cassert>
#include <iomanip>
#include "davidson.hpp"

template<typename T>
const size_t Davidson<T>::nsubreserve = 100;

template<typename T>
const dreal Davidson<T>::tol = 1e-6;

template<typename T>
void Davidson<T>::normalize(Vector<T> &vec)
{
  dreal norm = vec.norm();
  
  assert(norm > 0.0);
  
  vec = vec * static_cast<T>(1.0/norm);
}

template<typename T>
Vector<T> Davidson<T>::random_vector()
{
  size_t ndim = operate.get_ndim();
  
  Vector<T> vecrand(ndim);
  
  for(size_t idx = 0; idx < ndim; ++idx)
    vecrand[idx] = static_cast<T>( rand() / dreal(RAND_MAX) - 0.5 );
  
  normalize(vecrand);
  
  return vecrand;
}

template<typename T>
void Davidson<T>::orthogonalize(const std::vector< Vector<T> > &vecarr, Vector<T> &vec)
{
  Vector<T> veccopy = vec;
  
  for(size_t i = 0; i < vecarr.size(); ++i) {
    T coef = veccopy * vecarr[i];
    vec  = vec - vecarr[i] * coef;
  }
  
  normalize(vec);
}

template<typename T>
void Davidson<T>::expand_basis()
{
  if(basis.size() < 1) {
    basis.reserve(nsubreserve);
    cvec = random_vector();
  } else {
    orthogonalize(basis, cvec);
  }
  
  basis.push_back(cvec);
}

template<typename T>
void Davidson<T>::build_matred()
{
  size_t   nbasis = basis.size();
  SqMat<T> matsav = matred;
  
  matred.resize(nbasis);
  
  // Copy saved block
  for(size_t i = 0; i < matsav.size(); ++i)
    for(size_t j = 0; j < matsav.size(); ++j)
      matred.at(j, i) = matsav.at(j, i);
  // Augmented columns
  for(size_t i = matsav.size(); i < matred.size(); ++i) {
    Vector<T> oveci = basis[i];
    operate(oveci);
    for(size_t j = 0; j < matred.size(); ++j)
      matred.at(j, i) = basis[j] * oveci;
  }
  // Augmented rows
  for(size_t i = 0; i < matred.size(); ++i) {
    Vector<T> oveci = basis[i];
    operate(oveci);
    for(size_t j = matsav.size(); j < matred.size(); ++j)
      matred.at(j, i) = basis[j] * oveci;
  }
}

template<typename T>
void Davidson<T>::build_eigvec(const Vector<T> &vecred, Vector<T> &eigvec)
{
  assert(vecred.size() == basis.size());
  
  eigvec.clear();
  
  eigvec.resize(operate.get_ndim());
  
  for(size_t idx = 0; idx < vecred.size(); ++idx)
    eigvec = eigvec + basis[idx] * vecred[idx];
}

template<typename T>
void Davidson<T>::solve(dreal &eigval, Vector<T> &eigvec)
{
  SqMat<T>      mattmp;
  Vector<dreal> eig;
  Vector<T>     vecred, oeigvec;
  std::vector< Vector<T> > vecarr;
  
  do {
    
    expand_basis();
    
    build_matred();
    
    mattmp = matred;
    
    SqMat<T>::diagonalize(0, 0, mattmp, eig);
    Matrix<T>::vectorize(0, 0, mattmp, vecarr);
    
    eigval = eig[0];
    build_eigvec(vecarr[0], eigvec);
    
    oeigvec = eigvec;
    operate(oeigvec);
    
    rvec = oeigvec - eigvec * static_cast<T>(eigval);
    
    cvec = rvec;
    precond.set_rho(eigval);
    precond(cvec);
    
    std::cout << std::scientific << basis.size() << " " << rvec.norm() << std::endl;
    
  } while(rvec.norm() > tol && basis.size() < operate.get_ndim());
  
  std::cout << "Final size: " << basis.size() << std::endl;
} 

template class Davidson<dreal>;
template class Davidson<dcmplx>;
