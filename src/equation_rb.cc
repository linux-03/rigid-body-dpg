#ifndef EQUATION_RB_CC
#define EQUATION_RB_CC

#include "equation_rb.h"
#include "Newton.h"
#include "autodiff.h"
#include "chunked_jacobian.h"

namespace ASC_ode {

  RigidBodyEquation::RigidBodyEquation(RigidBody& rb, RigidBodySystem& rbs, double step_size): rb_(rb), rbs_(rbs), h_(step_size) {};

  size_t RigidBodyEquation::DimX() const {
    return this->dimX_;
  }

  size_t RigidBodyEquation::DimF() const {
    return this->dimF_;
  }

  template<typename T>
  Vector<T> RigidBodyEquation::FuncConstraint(const VectorView<T> x, bool old) const
  {
    Vector<T> res(12);
    res.setConstant(0);
    size_t constraint_count = 0;
    if (old)
    {
      Matrix<T> constr = rb_.Constraints();
      for (size_t i = 0; i < rb_.Beams().size(); i++)
      {
        Beam& bm = rbs_.Beams(rb_.Beams()[i]);
        res +=  constr.Rows(constraint_count, bm.NumberOfConstraints()/2).transpose() * x.segment(bm.StartIndex(), bm.NumberOfConstraints()/2);
        constraint_count += bm.NumberOfConstraints()/2;
      }
    } else {
      Matrix<T> constr = rbs_.JacobianConstraint(x, rb_.Index());
      
      for (size_t i = 0; i < rb_.Beams().size(); i++)
      {
        Beam& bm = rbs_.Beams(rb_.Beams()[i]);
        res += constr.Rows(constraint_count, bm.NumberOfConstraints()/2).transpose() * x.segment(bm.StartIndex() + bm.NumberOfConstraints()/2, bm.NumberOfConstraints()/2);
        constraint_count += bm.NumberOfConstraints()/2;
      }
    }

    return res;
  }


  template<typename T>
  void RigidBodyEquation::Residual(VectorView<T> x, VectorView<T> f) const
  {
    Vector<T> force_old = rb_.Force();
    Vector<T> vel_con_old = FuncConstraint(x, true);

    Vector<T> force_new = rbs_.PotentialGradient(x, rb_.Index());
    Vector<T> vel_con_new = FuncConstraint(x, false);

    // q_trans 3 q_rot 9 p_half 6 p 6 v 6 lambda 1 mu 1

    Vector<T> p = x.segment(dim_per_body()*rb_.Index() + 24, 6);
    Vector<T> v = x.segment(dim_per_body()*rb_.Index() + 12, 6);
    Vector<T> p_hat_new = x.segment(dim_per_body()*rb_.Index() + 18, 6);
    Vector<T> q_new = x.segment(dim_per_body()*rb_.Index(), 12);

    Matrix<T> Rmean = 0.5*(ToMatrix(q_new.segment(3, 9)) + rb_.q());
    Matrix<T> Rnew = ToMatrix(q_new.segment(3, 9));

    f.setConstant(0);

    f.segment(0, 3) = rb_.Mass() * v.segment(0, 3) - p.segment(0, 3);
    f.segment(3, 3) = v.segment(3, 3) - p.segment(3, 3);
    
    f.segment(6, 3) = q_new.segment(0, 3) - rb_.q_trans() - h_*v.segment(0, 3); // q_trans_n - q_trans - h_*p
    f.segment(9, 3) = skewSymmetricToVector(rb_.q().transpose() * ((ToMatrix(q_new.segment(3, 9)) - rb_.q()) - h_*Rmean*vectorToSkewSymmetric(v.segment(3, 3))));

    f.segment(12, 3) = p.segment(0, 3) - rb_.p_trans() - (h_/2)*(force_old.segment(0, 3) + vel_con_old.segment(0, 3));
    f.segment(15, 3) = skewSymmetricToVector( rb_.q().transpose() * ( Rmean * vectorToSkewSymmetric(v.segment(3, 3)) * rb_.TensorD() - rb_.q()* vectorToSkewSymmetric(rb_.p_skew()) - (h_/2)*(ToMatrix(force_old.segment(3, 9)) + ToMatrix(vel_con_old.segment(3, 9)))));

    f.segment(18, 3) = p_hat_new.segment(0, 3) - p.segment(0, 3) - (h_/2)*(force_new.segment(0, 3) + vel_con_new.segment(0, 3));
    f.segment(21, 3) = skewSymmetricToVector( Rnew.transpose() * (Rnew*vectorToSkewSymmetric(p_hat_new.segment(3, 3)) - Rmean * vectorToSkewSymmetric(v.segment(3, 3)) * rb_.TensorD() - (h_/2)*(ToMatrix(force_new.segment(3, 9)) + ToMatrix(vel_con_new.segment(3, 9)))));

    Matrix<T> eye(3, 3);
    eye.setConstant(0);
    eye(0, 0) = 1;
    eye(1, 1) = 1;
    eye(2, 2) = 1;
    Matrix<T> c = Rnew.transpose()*Rnew - eye;

    f(24) = c(0, 0);

    f(25) = c(1, 0);

    f(26) = c(2, 0);
    
    f(27) = c(1, 1);
    
    f(28) = c(2, 1);

    f(29) = c(2, 2);
  }

  void RigidBodyEquation::Evaluate(VectorView<double> x, VectorView<double> f) const
  {
    Residual<double>(x, f);
  }

  // The residual of body b reads
  //   - all 30 unknowns of b,
  //   - the 12 q of every free body that shares a beam or a spring with b,
  //   - the multipliers of the beams attached to b.
  // Every other column of its Jacobian rows is zero, so only these are seeded.
  void RigidBodyEquation::EvaluateDeriv(VectorView<double> x, MatrixView<double> df) const
  {
    df.setConstant(0);

    ColumnSet cols;
    cols.AddRange(dim_per_body()*rb_.Index(), dim_per_body());
    auto addQ = [&](Connector& c) {
      if (c.Type() != ConnectorType::FIX) cols.AddRange(dim_per_body()*c.BodyIndex(), 12);
    };
    for (size_t i = 0; i < rb_.Beams().size(); i++) {
      Beam& bm = rbs_.Beams(rb_.Beams()[i]);
      addQ(bm.ConnectorA());
      addQ(bm.ConnectorB());
      cols.AddRange(bm.StartIndex(), bm.NumberOfConstraints());
    }
    for (size_t k = 0; k < rbs_.NumSprings(); k++) {
      Spring& spr = rbs_.Springs(k);
      auto on_this = [&](Connector& c) { return c.Type() != ConnectorType::FIX && c.BodyIndex() == rb_.Index(); };
      if (on_this(spr.ConnectorA()) || on_this(spr.ConnectorB())) {
        addQ(spr.ConnectorA());
        addQ(spr.ConnectorB());
      }
    }

    ChunkedJacobian(x, cols.Finish(), dimF_,
                    [this](auto xd, auto fd) { Residual(xd, fd); }, df);
  }
}

#endif