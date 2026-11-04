#ifndef NONLINFUNC_H
#define NONLINFUNC_H

#include <mutex>
#include <thread>
#include <memory>
#include <vector.h>
#include <matrix.h>

using namespace ASC_bla;
using namespace std;



namespace ASC_ode
{


  class NonlinearFunction
  {
  public:
    virtual ~NonlinearFunction() = default;
    virtual size_t DimX() const = 0;
    virtual size_t DimF() const = 0;
    virtual void Evaluate (VectorView<double> x, VectorView<double> f) const = 0;
    virtual void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const = 0;
  };


  class IdentityFunction : public NonlinearFunction
  {
    size_t n;
  public:
    IdentityFunction (size_t _n);
    size_t DimX() const override;
    size_t DimF() const override;
    void Evaluate (VectorView<double> x, VectorView<double> f) const override;
    
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override;
  };

  class StackedFunction_large_input : public NonlinearFunction
  {
    vector<shared_ptr<NonlinearFunction>> _functions;
    public:
    StackedFunction_large_input();
    void addFunction(shared_ptr<NonlinearFunction> func);

    size_t DimX() const override;

    size_t DimF() const override;

    void Evaluate (VectorView<double> x, VectorView<double> f) const override;

    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override;
  };

  class StackedFunction : public NonlinearFunction
  {
    std::vector<std::shared_ptr<NonlinearFunction>> _functions;
    public:
    StackedFunction();
    void addFunction(std::shared_ptr<NonlinearFunction> func);
    size_t DimX() const override;

    size_t DimF() const override;

    void Evaluate (VectorView<double> x, VectorView<double> f) const override;

    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override;
  };

  void dNumeric(const NonlinearFunction& f, VectorView<double> x, MatrixView<double> df);
}
#endif