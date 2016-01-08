#ifndef DAVIDSON_HPP
#define DAVIDSON_HPP

#include <vector>
#include "precision.hpp"
#include "vector.hpp"
#include "sqmat.hpp"
#include "operator.hpp"

template<typename T>
class Davidson {
  
public:
  
  Davidson(Operator<T> &opr, Preconditioner<T> &pcd):
    operate(opr), precond(pcd), rvec(), cvec(), matred(), basis()
  { assert(opr.get_ndim() == pcd.get_ndim()); }
  
  void solve(dreal &eigval, Vector<T> &eigvec);
  
private:
  
  static const size_t nsubreserve;
  static const dreal  tol;
  
  Operator<T>       &operate;
  Preconditioner<T> &precond;
  
  Vector<T> rvec;
  Vector<T> cvec;
  SqMat<T>  matred;
  std::vector< Vector<T> > basis;
  
  Vector<T> random_vector();
  
  void normalize(Vector<T> &vec);
  
  void orthogonalize(const std::vector< Vector<T> > &vecarr, Vector<T> &vec);
  
  void expand_basis();
  
  void build_matred();
  
  void build_eigvec(const Vector<T> &vecred, Vector<T> &eigvec);
};

#endif
