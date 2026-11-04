#ifndef Newton_cc
#define Newton_cc

#include "Newton.h"
#include "iostream"

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ASC_ode {

  namespace {

    // Dense LU factorization with partial pivoting, stored row by row.
    class LU
    {
      size_t n_ = 0;
      std::vector<double> a_;
      std::vector<size_t> perm_;

    public:
      void factor(MatrixView<double> m)
      {
        n_ = m.Height();
        a_.resize(n_*n_);
        perm_.resize(n_);
        for (size_t i = 0; i < n_; i++)
        {
          perm_[i] = i;
          for (size_t j = 0; j < n_; j++)
            a_[i*n_ + j] = m(i, j);
        }

        double* a = a_.data();
        for (size_t k = 0; k < n_; k++)
        {
          size_t pivot = k;
          for (size_t i = k + 1; i < n_; i++)
            if (std::fabs(a[i*n_ + k]) > std::fabs(a[pivot*n_ + k]))
              pivot = i;
          if (a[pivot*n_ + k] == 0.0)
            throw std::domain_error("Matrix singular");
          if (pivot != k)
          {
            std::swap(perm_[pivot], perm_[k]);
            for (size_t j = 0; j < n_; j++)
              std::swap(a[pivot*n_ + j], a[k*n_ + j]);
          }
          const double* row_k = a + k*n_;
          for (size_t i = k + 1; i < n_; i++)
          {
            double* row_i = a + i*n_;
            if (row_i[k] == 0.0)
              continue;
            const double factor = row_i[k] / row_k[k];
            row_i[k] = factor;
            for (size_t j = k + 1; j < n_; j++)
              row_i[j] -= factor*row_k[j];
          }
        }
      }

      // Solves A*out = b.
      void solve(VectorView<double> b, std::vector<double>& out) const
      {
        const double* a = a_.data();
        out.resize(n_);
        for (size_t i = 0; i < n_; i++)
        {
          double sum = b(perm_[i]);
          for (size_t j = 0; j < i; j++)
            sum -= a[i*n_ + j]*out[j];
          out[i] = sum;
        }
        for (size_t i = n_; i-- > 0;)
        {
          double sum = out[i];
          for (size_t j = i + 1; j < n_; j++)
            sum -= a[i*n_ + j]*out[j];
          out[i] = sum / a[i*n_ + i];
        }
      }
    };

  }

  // Full Newton: the Jacobian is recomputed at every iteration. Reusing it would make
  // each step cheaper, but the leftover error of each solve accumulates, and energy and
  // momentum then drift far more over long runs. The tolerance is relative to the size
  // of the state, so it works for light and heavy bodies alike.
  void NewtonSolver (NonlinearFunction& func, VectorView<double> x, double tol, int maxsteps, std::function<void(int,double,VectorView<double>)> callback_)
  {
    Vector<double> res(func.DimF());
    res.setConstant(0);
    Matrix<double> fprime(func.DimF(), func.DimX());
    fprime.setConstant(0);
    LU lu;
    std::vector<double> update;

    for (int i = 0; i < maxsteps; i++)
      {
        func.Evaluate(x, res);
        const double err = res.norm();

        func.EvaluateDeriv(x, fprime);
        lu.factor(fprime);
        lu.solve(res, update);
        for (size_t k = 0; k < update.size(); k++)
          x(k) -= update[k];

        if (callback_)
          callback_(i, err, x);
        if (err < tol*(1.0 + x.norm())) return;
      }
    throw std::domain_error("Newton did not converge");
  }
}

#endif
