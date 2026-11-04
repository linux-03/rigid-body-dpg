#ifndef EQUATION_RB_CC
#define EQUATION_RB_CC

#include "equation_rb.h"
#include "Newton.h"
#include "autodiff.h"

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


  void RigidBodyEquation::Evaluate(VectorView<double> x, VectorView<double> f) const
  {
    f.setConstant(0);

    Vector<double> force_old = rb_.Force();
    Vector<double> vel_con_old = FuncConstraint(x, true);

    Vector<double> force_new = rbs_.PotentialGradient(x, rb_.Index());
    Vector<double> vel_con_new = FuncConstraint(x, false);

    // q_trans 3 q_rot 9 p_half 6 p 6 v 6 lambda 1 mu 1

    Vector<double> p = x.segment(dim_per_body()*rb_.Index() + 24, 6);
    Vector<double> v = x.segment(dim_per_body()*rb_.Index() + 12, 6);
    Vector<double> p_hat_new = x.segment(dim_per_body()*rb_.Index() + 18, 6);
    Vector<double> q_new = x.segment(dim_per_body()*rb_.Index(), 12);

    Matrix<double> Rmean = 0.5*(ToMatrix(q_new.segment(3, 9)) + rb_.q());
    Matrix<double> Rnew = ToMatrix(q_new.segment(3, 9));

    f.setConstant(0);

    f.segment(0, 3) = rb_.Mass() * v.segment(0, 3) - p.segment(0, 3);
    f.segment(3, 3) = v.segment(3, 3) - p.segment(3, 3);
    
    f.segment(6, 3) = q_new.segment(0, 3) - rb_.q_trans() - h_*v.segment(0, 3); // q_trans_n - q_trans - h_*p
    f.segment(9, 3) = skewSymmetricToVector(rb_.q().transpose() * ((ToMatrix(q_new.segment(3, 9)) - rb_.q()) - h_*Rmean*vectorToSkewSymmetric(v.segment(3, 3))));

    f.segment(12, 3) = p.segment(0, 3) - rb_.p_trans() - (h_/2)*(force_old.segment(0, 3) + vel_con_old.segment(0, 3));
    f.segment(15, 3) = skewSymmetricToVector( rb_.q().transpose() * ( Rmean * vectorToSkewSymmetric(v.segment(3, 3)) * rb_.TensorD() - rb_.q()* vectorToSkewSymmetric(rb_.p_skew()) - (h_/2)*(ToMatrix(force_old.segment(3, 9)) + ToMatrix(vel_con_old.segment(3, 9)))));

    f.segment(18, 3) = p_hat_new.segment(0, 3) - p.segment(0, 3) - (h_/2)*(force_new.segment(0, 3) + vel_con_new.segment(0, 3));
    f.segment(21, 3) = skewSymmetricToVector( Rnew.transpose() * (Rnew*vectorToSkewSymmetric(p_hat_new.segment(3, 3)) - Rmean * vectorToSkewSymmetric(v.segment(3, 3)) * rb_.TensorD() - (h_/2)*(ToMatrix(force_new.segment(3, 9)) + ToMatrix(vel_con_new.segment(3, 9)))));

    Matrix<double> eye(3, 3);
    eye.setConstant(0);
    eye(0, 0) = 1;
    eye(1, 1) = 1;
    eye(2, 2) = 1;
    Matrix<double> c = Rnew.transpose()*Rnew - eye;

    f(24) = c(0, 0);

    f(25) = c(1, 0);

    f(26) = c(2, 0);
    
    f(27) = c(1, 1);
    
    f(28) = c(2, 1);

    f(29) = c(2, 2);
  }

  void RigidBodyEquation::EvaluateDeriv(VectorView<double> x, MatrixView<double> df)  const {

    df.setConstant(0);

    Vector<AutoDiff<30, double>> x_diff = x;
    size_t prev_index = 0;
    Vector<AutoDiff<30, double>> f_diff(30);

    for (size_t bd_index = 0; bd_index < rbs_.NumBodies(); bd_index++) {

      for (size_t j = 0; j < 30; j++) {
        x_diff(dim_per_body()*prev_index + j).DValue(j) = 0;
        x_diff(dim_per_body()*bd_index + j).DValue(j) = 1;
      }
      prev_index = bd_index;

      Vector<AutoDiff<30, double>> force_old = rb_.Force();
      Vector<AutoDiff<30, double>> vel_con_old = FuncConstraint(x_diff, true);
      

      Vector<AutoDiff<30, double>> force_new = rbs_.PotentialGradient(x_diff, rb_.Index());
      Vector<AutoDiff<30, double>> vel_con_new = FuncConstraint(x_diff, false);

      Vector<AutoDiff<30, double>> p = x_diff.segment(dim_per_body()*rb_.Index() + 24, 6);
      Vector<AutoDiff<30, double>> v = x_diff.segment(dim_per_body()*rb_.Index() + 12, 6);
      Vector<AutoDiff<30, double>> p_hat_new = x_diff.segment(dim_per_body()*rb_.Index() + 18, 6);
      Vector<AutoDiff<30, double>> q_new = x_diff.segment(dim_per_body()*rb_.Index(), 12);

      Matrix<AutoDiff<30, double>> Rmean = 0.5*(ToMatrix(q_new.segment(3, 9)) + rb_.q());
      Matrix<AutoDiff<30, double>> Rnew = ToMatrix(q_new.segment(3, 9));

      // q_trans 3 q_rot 9 p_half 6 p 6 v 6 lambda 1 mu 1
      
      
      f_diff.segment(0, 3) = rb_.Mass() * v.segment(0, 3) - p.segment(0, 3);
      f_diff.segment(3, 3) = v.segment(3, 3) - p.segment(3, 3);
    
      f_diff.segment(6, 3) = q_new.segment(0, 3) - rb_.q_trans() - h_*v.segment(0, 3); // q_trans_n - q_trans - h_*p
      f_diff.segment(9, 3) = skewSymmetricToVector(rb_.q().transpose() * ((ToMatrix(q_new.segment(3, 9)) - rb_.q()) - h_*Rmean*vectorToSkewSymmetric(v.segment(3, 3))));

      f_diff.segment(12, 3) = p.segment(0, 3) - rb_.p_trans() - (h_/2)*(force_old.segment(0, 3) + vel_con_old.segment(0, 3));  
      f_diff.segment(15, 3) = skewSymmetricToVector( rb_.q().transpose() * ( Rmean * vectorToSkewSymmetric(v.segment(3, 3)) * rb_.TensorD() - rb_.q()* vectorToSkewSymmetric(rb_.p_skew()) - (h_/2)*(ToMatrix(force_old.segment(3, 9)) + ToMatrix(vel_con_old.segment(3, 9)))));

      f_diff.segment(18, 3) = p_hat_new.segment(0, 3) - p.segment(0, 3) - (h_/2)*(force_new.segment(0, 3) + vel_con_new.segment(0, 3));
      f_diff.segment(21, 3) = skewSymmetricToVector( Rnew.transpose() * (Rnew*vectorToSkewSymmetric(p_hat_new.segment(3, 3)) - Rmean * vectorToSkewSymmetric(v.segment(3, 3)) * rb_.TensorD() - (h_/2)*(ToMatrix(force_new.segment(3, 9)) + ToMatrix(vel_con_new.segment(3, 9)))));


      Matrix<AutoDiff<30, double>> eye(3, 3);
      eye.setConstant(0);
      eye(0, 0) = 1;
      eye(1, 1) = 1;
      eye(2, 2) = 1;
      Matrix<AutoDiff<30, double>> c = Rnew.transpose()*Rnew - eye;

      f_diff(24) = c(0, 0);

      f_diff(25) = c(1, 0);

      f_diff(26) = c(2, 0);
      
      f_diff(27) = c(1, 1);
      
      f_diff(28) = c(2, 1);

      f_diff(29) = c(2, 2);

      for (size_t j = 0; j < 30; j++) {
        df.Row(j).segment(dim_per_body()*bd_index, dim_per_body()) = f_diff(j); 
      }
    }

    prev_index = 0;

    Vector<AutoDiff<2, double>> x_bm_diff = x;
    Vector<AutoDiff<2, double>> f_bm_diff(30);

    for (size_t bm_index = 0; bm_index < rbs_.NumConstraints()/2; bm_index++) {

      x_bm_diff(dim_per_body()*rbs_.NumBodies() + 2*prev_index).DValue(0) = 0;
      x_bm_diff(dim_per_body()*rbs_.NumBodies() + 2*prev_index + 1).DValue(1) = 0;

      x_bm_diff(dim_per_body()*rbs_.NumBodies() + 2*bm_index).DValue(0) = 1;
      x_bm_diff(dim_per_body()*rbs_.NumBodies() + 2*bm_index + 1).DValue(1) = 1;
      prev_index = bm_index;

      Vector<AutoDiff<2, double>> force_old = rb_.Force();
      Vector<AutoDiff<2, double>> vel_con_old = FuncConstraint(x_bm_diff, true);

      Vector<AutoDiff<2, double>> force_new = rbs_.PotentialGradient(x_bm_diff, rb_.Index());
      Vector<AutoDiff<2, double>> vel_con_new = FuncConstraint(x_bm_diff, false);

      Vector<AutoDiff<2, double>> p = x_bm_diff.segment(dim_per_body()*rb_.Index() + 24, 6);
      Vector<AutoDiff<2, double>> v = x_bm_diff.segment(dim_per_body()*rb_.Index() + 12, 6);
      Vector<AutoDiff<2, double>> p_hat_new = x_bm_diff.segment(dim_per_body()*rb_.Index() + 18, 6);
      Vector<AutoDiff<2, double>> q_new = x_bm_diff.segment(dim_per_body()*rb_.Index(), 12);

      Matrix<AutoDiff<2, double>> Rmean = 0.5*(ToMatrix(q_new.segment(3, 9)) + rb_.q());
      Matrix<AutoDiff<2, double>> Rnew = ToMatrix(q_new.segment(3, 9));

      
      
      f_bm_diff.segment(0, 3) = rb_.Mass() * v.segment(0, 3) - p.segment(0, 3);
      f_bm_diff.segment(3, 3) = v.segment(3, 3) - p.segment(3, 3);
    
      f_bm_diff.segment(6, 3) = q_new.segment(0, 3) - rb_.q_trans() - h_*v.segment(0, 3); // q_trans_n - q_trans - h_*p
      f_bm_diff.segment(9, 3) = skewSymmetricToVector(rb_.q().transpose() * ((ToMatrix(q_new.segment(3, 9)) - rb_.q()) - h_*Rmean*vectorToSkewSymmetric(v.segment(3, 3))));

      f_bm_diff.segment(12, 3) = p.segment(0, 3) - rb_.p_trans() - (h_/2)*(force_old.segment(0, 3) + vel_con_old.segment(0, 3));
      f_bm_diff.segment(15, 3) = skewSymmetricToVector( rb_.q().transpose() * ( Rmean * vectorToSkewSymmetric(v.segment(3, 3)) * rb_.TensorD() - rb_.q()* vectorToSkewSymmetric(rb_.p_skew()) - (h_/2)*(ToMatrix(force_old.segment(3, 9)) + ToMatrix(vel_con_old.segment(3, 9)))));

      f_bm_diff.segment(18, 3) = p_hat_new.segment(0, 3) - p.segment(0, 3) - (h_/2)*(force_new.segment(0, 3) + vel_con_new.segment(0, 3));
      f_bm_diff.segment(21, 3) = skewSymmetricToVector( Rnew.transpose() * (Rnew*vectorToSkewSymmetric(p_hat_new.segment(3, 3)) - Rmean * vectorToSkewSymmetric(v.segment(3, 3)) * rb_.TensorD() - (h_/2)*(ToMatrix(force_new.segment(3, 9)) + ToMatrix(vel_con_new.segment(3, 9)))));

      for (size_t j = 0; j < 30; j++) {
        df.Row(j).segment(dim_per_body()*rbs_.NumBodies() + bm_index*2, 2) = f_bm_diff(j); 
      }
    }    
  }
}

#endif