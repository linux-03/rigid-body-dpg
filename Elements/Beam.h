#ifndef BEAM_H
#define BEAM_H

#include <vector.h>
#include "Connector.h"

using namespace ASC_bla;

class Beam {
    size_t index_;
    double length_;
    double lambda_;
    double mu_;
    size_t start_index_ = 0;
    size_t number_of_constraints_ = 2;
    std::pair<Vector<double>, Vector<double>> constraint_plane_b_ = {Vector<double>(3), Vector<double>(3)};
    std::pair<Vector<double>, Vector<double>> constraint_plane_a_ = {Vector<double>(3), Vector<double>(3)};
    Vector<double> axis_b_ = {0, 0, 0};
    Vector<double> axis_a_ = {0, 0, 0};
    Vector<double> lambda_vec_;
    Vector<double> mu_vec_;
    Connector connector_a_;
    Connector connector_b_;

  public:
    Beam (Connector connector_a, Connector connector_b);
    Beam (const Beam& other);

    size_t& Index();
    double& Length();

    size_t BodyIndexA();
    size_t BodyIndexB();

    size_t& StartIndex();

    Vector<double>& LambdaVec();
    Vector<double>& MuVec();

    Connector& ConnectorA();
    Connector& ConnectorB();

    template<typename T>
    Vector<T> BodyFramePositionA(const VectorView<T> q) {
      return this->connector_a_.BodyFramePosition(q);
    }

    template<typename T>
    Vector<T> BodyFramePositionB(const VectorView<T> q) {
      return this->connector_b_.BodyFramePosition(q);
    }

    template<typename T>
    Vector<T> AbsPositionA(const VectorView<T> q) {
      return this->connector_a_.AbsPosition(q);
    }

    template<typename T>
    Vector<T> AbsPositionB(const VectorView<T> q) {
      return this->connector_b_.AbsPosition(q);
    }

    template<typename T>
    Vector<T> Constraint(const VectorView<T> q_a, const VectorView<T> q_b) {
      Vector<T> res(this->number_of_constraints_/2);
      //if( !(this->connector_a_.Type() == ConnectorType::SPHERICAL) && !(this->connector_b_.Type() == ConnectorType::SPHERICAL)) {
      Vector<T> q = AbsPositionA(q_a) - AbsPositionB(q_b);
      
      res(0) = 0.5*(q.squaredNorm() - this->length_*this->length_);
        
      //}
      size_t constraint_count = 1;
      if (this->connector_a_.Type() == ConnectorType::SPHERICAL) {
        Matrix<T> R = ToMatrix(q_a.segment(3, 9));

        Vector<T> test_1 = R * this->axis_a_;
        Vector<T> test_2 = R * this->constraint_plane_a_.second;
        Vector<T> test_3 = R * this->constraint_plane_a_.first;
        
        Vector<T> test = ApplyTransformation(q_a, this->axis_a_) - AbsPositionB(q_b);
        
        res(constraint_count) = test * test_2;
        res(constraint_count + 1) = test * test_3;

        constraint_count +=2;
        
      }
      if (this->connector_b_.Type() == ConnectorType::SPHERICAL) {
        Matrix<T> R = ToMatrix(q_b.segment(3, 9));

        Vector<T> test_1 = R * this->axis_b_;
        Vector<T> test_2 = R * this->constraint_plane_b_.second;
        Vector<T> test_3 = R * this->constraint_plane_b_.first;
        
        Vector<T> test = ApplyTransformation(q_b, this->axis_b_) - AbsPositionA(q_a);
        
        res(constraint_count) = test * test_2;
        res(constraint_count + 1) = test * test_3;
        
      }
      return res;
    }

    size_t NumberOfConstraints();
};

#endif