#ifndef OPERATOR_HPP
#define OPERATOR_HPP

template<typename T>
class Operator {
  
public:
  
  Operator(size_t n): ndim(n) {}
  
  virtual ~Operator() {}
  
  size_t get_ndim() const
  { return ndim; }
  
  virtual void operator()(Vector<T> &vec) = 0;
  
protected:
  const size_t ndim;
};

template<typename T>
class Preconditioner : public Operator<T> {
  
public:
  
  Preconditioner(size_t n): Operator<T>(n), rho(0) {}
  
  virtual ~Preconditioner() {}
  
  virtual void operator()(Vector<T> &vec) = 0;
  
  void set_rho(dreal x)
  { rho = x; }
  
protected:
  dreal rho;
};

#endif
