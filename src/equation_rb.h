#ifndef EQUATION_RB_H
#define EQUATION_RB_H

#include "nonlinfunc.h"
#include "rigid_body_system.h"

namespace ASC_ode {

  class RigidBodyEquation: public NonlinearFunction
  {
      size_t dimX_ = 30;
      size_t dimF_ = 30;
      RigidBody& rb_;
      RigidBodySystem& rbs_;
      double h_;

    public:
      RigidBodyEquation(RigidBody& rb, RigidBodySystem& rbs, double step_size);
      size_t DimX() const;
      size_t DimF() const;
      void Evaluate (VectorView<double> x, VectorView<double> f) const;
      void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const;

      // The body's 30 equations, for double values and for AutoDiff values alike.
      template<typename T>
      void Residual(VectorView<T> x, VectorView<T> f) const;
      
      template<typename T>
      Vector<T> FuncConstraint(const VectorView<T> x, bool old) const;
  };

  //void simulate(RigidBody& rb, double step_size, size_t steps_);

}

#endif