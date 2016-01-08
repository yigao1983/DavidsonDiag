#include <iostream>
#include "vector.hpp"
#include "matrix.hpp"
#include "sqmat.hpp"
#include "davidson.hpp"

using namespace std;

class Functor : public Operator<dcmplx> {
  
public:
  
  Functor(const SqMatDcmplx &sqmat):
    Operator<dcmplx>(sqmat.size()), omat(sqmat)
  {}
  
  virtual ~Functor() {}
  
  virtual void operator()(VectorDcmplx &vec)
  {
    assert(vec.size() == Operator<dcmplx>::ndim);
    vec = omat * vec;
  }
  
protected:
  const SqMatDcmplx &omat;
};

class PreFunctor : public Preconditioner<dcmplx> {
  
public:
  
  PreFunctor(const SqMatDcmplx &sqmat):
    Preconditioner<dcmplx>(sqmat.size()), omat(sqmat)
  {}
  
  virtual ~PreFunctor() {}
  
  virtual void operator()(VectorDcmplx &vec)
  {
    assert(vec.size() == Operator<dcmplx>::ndim);
    
    for(size_t idx = 0; idx < vec.size(); ++idx) {
      dcmplx diff = rho - omat.at(idx, idx);
      if(abs(diff) > 0.0) vec[idx] = vec[idx] / diff;
    }
  }
  
protected:
  const SqMatDcmplx &omat;
};

int main()
{
  const size_t n = 1000;
  SqMatDcmplx sqmat(n);
  
  for(size_t i = 0; i < sqmat.size(); ++i)
    for(size_t j = 0; j <= i; ++j) {
      sqmat.at(j, i) = static_cast<dcmplx>(rand() / dreal(RAND_MAX) - 0.5);
      sqmat.at(i, j) = sqmat.at(j, i);
    }
  
  Functor opr(sqmat);
  PreFunctor pcd(sqmat);
  Davidson<dcmplx> davidson(opr, pcd);
  
  dreal x;
  Vector<dcmplx> vec;
  davidson.solve(x, vec);
  
  cout << "Davidson: x " << x << endl;
  
  VectorDreal eig;
  
  SqMatDcmplx::diagonalize(0, 0, sqmat, eig);
  
  vector<VectorDcmplx> vecarr;
  MatDcmplx::vectorize(0, 0, sqmat, vecarr);
  
  cout << "Exact: x " << eig[0] << endl;
  
  cout << "accuracy: " << (vecarr[0]-vec).norm() << endl;
  
  return 0;
}
