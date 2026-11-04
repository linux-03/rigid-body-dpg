#ifndef RIGID_BODY_SYSTEM_H
#define RIGID_BODY_SYSTEM_H

#include "rigid_body.h"
#include "Beam.h"
#include "Spring.h"
#include "Connector.h"

using namespace std;

class RigidBodySystem {
    size_t num_beams_;
    size_t num_springs_;
    size_t num_constraints_ = 0;
    std::vector<std::reference_wrapper<RigidBody>> bodies_;
    std::vector<std::reference_wrapper<Beam>> beams_;
    std::vector<std::reference_wrapper<Spring>> springs_;
    Vector<double> gravity_ = {0, -9.81, 0};

  public:
    RigidBodySystem();

    template<typename T>
    T Potential(const VectorView<T> x);

    /**
    * @brief Computes the applied force on one body, which is minus the gradient of the potential
    *
    * @param x The full position vector containing each position and momentum of each body
    * @param body_index The index of the body to get the gradient fo
    * @return 12 Dimensional force
    */
    template<typename T>
    Vector<T> PotentialGradient(const VectorView<T> x, size_t body_index)
    {
      Vector<T> pot_grad(12);
      pot_grad.setConstant(0);

      pot_grad.segment(0, 3) = this->Bodies(body_index).Mass()*gravity_;

      for (size_t k = 0; k < springs_.size(); k++)
      {
        Spring& spr = springs_[k];
        Connector& a = spr.ConnectorA();
        Connector& b = spr.ConnectorB();
        Vector<T> force = spr.ForceOnA(x.segment(a.BodyIndex()*dim_per_body(), 12), x.segment(b.BodyIndex()*dim_per_body(), 12));
        // Generalized force on (q, Q) from a force f at the body point s: f, and f_i*s_j.
        auto apply = [&](Connector& c, double sign) {
          if (c.Type() == ConnectorType::FIX || c.BodyIndex() != body_index) return;
          Vector<double> s = c.RefPosition();
          for (size_t i = 0; i < 3; i++)
          {
            pot_grad(i) = pot_grad(i) + sign*force(i);
            for (size_t j = 0; j < 3; j++)
              pot_grad(3 + 3*i + j) = pot_grad(3 + 3*i + j) + (sign*s(j))*force(i);
          }
        };
        apply(a, 1.0);
        apply(b, -1.0);
      }

      return pot_grad;
    }


    /**
    * @brief Computes the length constraint of a beam
    *
    * @param x The full position vector containing each position and momentum of each body
    * @param beam_index The index of the beam to have the constraint calculated
    * @return Value T of the length difference
    */
    template<typename T>
    Vector<T> Constraint(const VectorView<T> x, size_t beam_index)
    {
      Beam& bm = this->beams_[beam_index];
      return bm.Constraint(x.segment(bm.BodyIndexA()*dim_per_body(), 12), x.segment(bm.BodyIndexB()*dim_per_body(), 12));
    }


    /**
    * @brief Computes the Jacobien of the function of all constraint derived after a specific body
    *
    * @param x The full position vector containing each position and momentum of each body
    * @param body_index The index of the body to derive after
    * @return Jacobina Matrix with heigth according to the number of beams and width 12
    */
    template<typename T>
    Matrix<T> JacobianConstraint(const VectorView<T> x, size_t body_index)
    {
      Matrix<T> G(this->Bodies(body_index).ConstraintNumberHalf(), 12);
      G.setConstant(0);
      double constraint_counter = 0;
      for (size_t i = 0; i < this->Bodies(body_index).Beams().size(); i++)
      { 
        Beam bm = Beams(this->Bodies(body_index).Beams()[i]);
        size_t index_a = bm.BodyIndexA();
        size_t index_b = bm.BodyIndexB();
        G.Block(constraint_counter, 0, bm.NumberOfConstraints()/2, 12) = G_single_body_beam_n(bm, body_index, x.segment(index_a*dim_per_body(), 12), x.segment(index_b*dim_per_body(), 12));
        constraint_counter += bm.NumberOfConstraints()/2;
      }

      return G;
    }


    /**
    * @brief Computes the derivative of a beam length constrain with respect to a specific body.
    *
    * @param bm The Beam of which you want to have the derivative of the constraint
    * @param bd_index The index of the body to derive after
    * @param q_a The 12-dimesnional vector representing the body position of body A of the beam
    * @param q_b The 12-dimesnional vector representing the body position of body B of
    * @return 12-dimesional gradient
    */
    template<typename T>
    Matrix<T> G_single_body_beam_n(Beam& bm, size_t bd_index, VectorView<T> q_a, VectorView<T> q_b) {
    
      Matrix<T> res(bm.NumberOfConstraints()/2 , dim_per_motion());

      Vec<3, T> pos1 = bm.ConnectorA().RefPosition();
      Vec<3, T> pos2 = bm.ConnectorB().RefPosition();

      Vector<AutoDiff<12, T>> q_a_diff = q_a;
      Vector<AutoDiff<12, T>> q_b_diff = q_b;

      if ((bd_index == bm.BodyIndexA()) && !(bm.ConnectorA().Type() == ConnectorType::FIX)) {
        for(size_t i =0 ; i < 12; i++){
          q_a_diff(i).DValue(i) = 1;
        }
      } else {
        for(size_t i=0; i < 12; i++){
          q_b_diff(i).DValue(i) = 1;
        }
      }

      Vector<AutoDiff<12, T>> g_vec = bm.Constraint(q_a_diff, q_b_diff);
      for (size_t j = 0; j < res.Height(); j++) {
        for(size_t i = 0 ; i < 12; i++){
          res(j, i) = g_vec(j).DValue(i);
        }
      }
      return res;
    }


    /**
    * @brief Computes the VelocityConstraint for a specific beam
    *
    * @param x The full position vector containing each position and momentum of each body
    * @param beam_index The index of the beam to calculate for
    * @return T value of the constraint
    */
    template<typename T>
    Vector<T> Velocity_Constraint(const VectorView<T> x, size_t beam_index)
    { 
      Beam& bm = this->beams_[beam_index];
      Matrix<T> G_a(bm.NumberOfConstraints()/2 , dim_per_motion()); // Jacobian for one constraint and three translational coordinates
      Matrix<T> G_b(bm.NumberOfConstraints()/2 , dim_per_motion());

      G_a.setConstant(0);
      G_b.setConstant(0);

      size_t bd_index_a = bm.BodyIndexA();
      size_t bd_index_b = bm.BodyIndexB();

      G_a += G_single_body_beam_n(bm, bd_index_a, x.segment(bd_index_a*dim_per_body(), 12), x.segment(bd_index_b*dim_per_body(), 12));
      G_b += G_single_body_beam_n(bm, bd_index_b, x.segment(bd_index_a*dim_per_body(), 12), x.segment(bd_index_b*dim_per_body(), 12));
      
      // G * M^-1 * p for one body: velocity p/m, and Qdot = Q*hat(omega) with
      // omega = 2*J^-1*p_skew. Fixed connectors contribute no velocity.
      auto rate=[&](Matrix<T>& G, size_t bd){
        RigidBody& b=this->Bodies(bd);
        Vector<T> pt=x.segment(bd*dim_per_body()+18,3);
        Vector<T> phat=x.segment(bd*dim_per_body()+21,3);
        MatrixView<double> Jinv=b.InertiaInv();
        Vector<T> w(3); Vector<T> pt2(3);
        for(size_t k=0;k<3;k++){ pt2(k)=(1.0/b.Mass())*pt(k); T acc=0.0*phat(0); for(size_t l=0;l<3;l++) acc=acc+(2.0*Jinv(k,l))*phat(l); w(k)=acc; }
        Matrix<T> R=ToMatrix(x.segment(bd*dim_per_body()+3,9))*vectorToSkewSymmetric(VectorView<T>(w));
        Vector<T> res=G.Cols(0,3)*pt2 + G.Cols(3,9)*ToVector(MatrixView<T>(R));
        return res;
      };
      Vector<T> out(bm.NumberOfConstraints()/2); out.setConstant(0);
      if (bm.ConnectorA().Type()!=ConnectorType::FIX) { Vector<T> c=rate(G_a,bd_index_a); out=out+c; }
      if (bm.ConnectorB().Type()!=ConnectorType::FIX) { Vector<T> c=rate(G_b,bd_index_b); out=out+c; }
      return out;
    }

    Vector<double>& Gravity() { return gravity_; }
    size_t& NumBeams();
    size_t& NumSprings();
    size_t NumBodies() const;
    size_t BodyDimensions() const;
    size_t NumConstraints() const;
    
    RigidBody& Bodies(size_t i);
    Beam& Beams(size_t i);
    std::vector<Beam> Beams();
    std::vector<RigidBody> Bodies();
    Spring& Springs(size_t i);

    Vector<double> connectorPosition(Connector c);

    void SaveState(const VectorView<double> x);
    void ManageConstraints(const VectorView<double> x);
    Vector<double> ExpandState();
    

    void add(RigidBody& rb);
    void add(Beam& bm);
    void add(Spring& spr);
};

#endif