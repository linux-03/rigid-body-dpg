#ifndef NONLINFUNC_CC
#define NONLINFUNC_CC

#include "nonlinfunc.h"
#include <iostream>

using namespace std;



namespace ASC_ode
{
  //  return Map<MatrixXd>(m.block(rows_start, cols_start, rows_num, cols_num).data(), rows_num, cols_num);
  //}

  //  return Map<VectorXd>(v.data(), v.size());
  //}
  
  IdentityFunction::IdentityFunction (size_t _n) : n(_n) { } 
  size_t IdentityFunction::DimX() const { return n; }
  size_t IdentityFunction::DimF() const { return n; }
  void IdentityFunction::Evaluate (VectorView<double> x, VectorView<double> f) const {
    f = x;
  }
  
  void IdentityFunction::EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const {
    df.setConstant(0);
    df.diagonal(1);
  }
  
  StackedFunction_large_input::StackedFunction_large_input(){};
  void StackedFunction_large_input::addFunction(shared_ptr<NonlinearFunction> func)
  {
    _functions.push_back(func);
  }

  size_t StackedFunction_large_input::DimX() const {
    size_t sum = 0;
    for(auto ptr : _functions)
      sum += ptr->DimX();
    return sum;
  }

  size_t StackedFunction_large_input::DimF() const {
    size_t sum=0;
    for(auto ptr : _functions)
      sum += ptr->DimF();
    return sum;
  }

  void StackedFunction_large_input::Evaluate (VectorView<double> x, VectorView<double> f) const {
    size_t cursor_f=0;
    for(auto func : _functions){

      func->Evaluate(x, f.segment(cursor_f, func->DimF()));





      cursor_f+=func->DimF();
    }
  }

  void StackedFunction_large_input::EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const{
    size_t cursor_f=0;
    df.setConstant(0);
    for(size_t i=0; i <_functions.size();i++ ){
      auto func = _functions[i];
      //The Jacobian of the stacked function is a block diagonal matrix, consisting of the individual Jacobians
      func->EvaluateDeriv(x, df.Block(cursor_f, 0, func->DimF(), df.Width()));

      cursor_f+=func->DimF();
    }
  }
  
  StackedFunction::StackedFunction() {};
  void StackedFunction::addFunction(std::shared_ptr<NonlinearFunction> func){_functions.push_back(func);}
  size_t StackedFunction::DimX() const {
    size_t sum = 0;
    for(auto ptr : _functions)
      sum += ptr->DimX();
    return sum;
  }

  size_t StackedFunction::DimF() const {
    size_t sum=0;
    for(auto ptr : _functions)
      sum += ptr->DimF();
    return sum;
  }
  void StackedFunction::Evaluate (VectorView<double> x, VectorView<double> f) const{
    size_t cursor_x=0;
    size_t cursor_f=0;
    for(auto func : _functions){

      func->Evaluate(x.segment(cursor_x, func->DimX()), f.segment(cursor_f, func->DimF()));


      cursor_f+=func->DimF();
      cursor_x+=func->DimX();
    }
  }

  void StackedFunction::EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const{
    size_t cursor_x=0;
    size_t cursor_f=0;
    df.setConstant(0);
    for(int i=0;i<_functions.size();i++ ){
      auto func = _functions[i];
      //The Jacobian of the stacked function is a block diagonal matrix, consisting of the individual Jacobians

      func->EvaluateDeriv(x.segment(cursor_x, func->DimX()), df.Block(cursor_f, cursor_x, func->DimF(), func->DimX()));

      cursor_f+=func->DimF();
      cursor_x+=func->DimX();
    }
  }

  void dNumeric(const NonlinearFunction& f, VectorView<double> x, MatrixView<double> df){
    double eps = 1e-6;
    Vector<double> xl(x.Size()), xr(x.Size()), fl(f.DimF()), fr(f.DimF());

    for (size_t i = 0; i < x.Size(); i++)
      {
        fl.setConstant(0);
        fr.setConstant(0);
        xl = x;
        xl(i) -= eps;
        
        xr = x;
        xr(i) += eps;
        
        f.Evaluate(xl, fl);
        f.Evaluate(xr, fr);
        
        df.Col(i) = 1/(2*eps) * (fr-fl);

      }
  }
}
#endif