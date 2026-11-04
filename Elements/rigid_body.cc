#ifndef RIGID_BODY_CC
#define RIGID_BODY_CC

#include "rigid_body.h"

using namespace ASC_bla;

RigidBody::RigidBody(): force_(12), q_trans_(3), q_(3, 3), p_trans_(3),
p_skew_(3), p_half_trans_(3), p_half_skew_(3), v_trans_(3), v_skew_(3), mass_(1.0), inertia_(Diagonal(3, 1.0)), mass_matrix_(Diagonal(6, 1.0)), mass_matrix_inv_(6,6)
{
  q_trans_.setConstant(0);
  q_.setConstant(0);
  p_trans_.setConstant(0);
  p_skew_.setConstant(0);
  p_half_trans_.setConstant(0);
  p_half_skew_.setConstant(0);
  v_trans_.setConstant(0);
  v_skew_.setConstant(0);

  q_.diagonal(1);
  constraints_.reset(new Matrix<double>(0,12));
  recalcMassMatrix();
}

RigidBody::RigidBody(VectorView<double> q, VectorView<double> p): force_(12), q_trans_(3), q_(3, 3), p_trans_(3),
p_skew_(3), p_half_trans_(3), p_half_skew_(3), v_trans_(3), v_skew_(3), mass_(1.0), inertia_(Diagonal(3, 1.0)), mass_matrix_(Diagonal(6, 1.0)), mass_matrix_inv_(6,6)
{
  q_trans_ = q.segment(0, 3);
  q_ = ToMatrix(q.segment(3, 9));
  p_trans_ = p.segment(0, 3);
  p_skew_ = p.segment(3, 3);
  p_half_trans_.setConstant(0);
  p_half_skew_.setConstant(0);
  v_trans_.setConstant(0);
  v_skew_.setConstant(0);
  constraints_.reset(new Matrix<double>(0,12));
  recalcMassMatrix();
}

RigidBody::RigidBody(double mass, MatrixView<double> inertia): force_(12), q_trans_(3), q_(3, 3), p_trans_(3),
p_skew_(3), p_half_trans_(3), p_half_skew_(3), v_trans_(3), v_skew_(3), mass_(mass), inertia_(inertia), mass_matrix_(6,6), mass_matrix_inv_(6,6)
{
  q_trans_.setConstant(0);
  q_.setConstant(0);
  p_trans_.setConstant(0);
  p_skew_.setConstant(0);
  p_half_trans_.setConstant(0);
  p_half_skew_.setConstant(0);
  v_trans_.setConstant(0);
  v_skew_.setConstant(0);

  q_.diagonal(1);
  constraints_.reset(new Matrix<double>(0,12));
  recalcMassMatrix();
}


RigidBody::RigidBody(const RigidBody& other): force_(12), q_trans_(3), q_(3, 3), p_trans_(3),
p_skew_(3), p_half_trans_(3), p_half_skew_(3), v_trans_(3), v_skew_(3), mass_(other.mass_), inertia_(other.inertia_), mass_matrix_(other.mass_matrix_), mass_matrix_inv_(other.mass_matrix_inv_) {
  q_trans_ = other.q_trans_;
  q_ = other.q_;
  p_trans_ = other.p_trans_;
  p_skew_ = other.p_skew_;
  p_half_trans_ = other.p_half_trans_;
  p_half_skew_ = other.p_half_skew_;
  v_trans_ = other.v_trans_;
  v_skew_ = other.v_skew_;
  lambda_ = other.lambda_;
  mu_ = other.mu_;
  beams_ = other.beams_;
  body_index_ = other.body_index_;
  constraint_number_half_ = other.constraint_number_half_;
  if (other.constraints_) constraints_.reset(new Matrix<double>(*other.constraints_));
  force_ = other.force_;
  inertia_inv_ = other.inertia_inv_;
  vertices_ = other.vertices_;
  normals_ = other.normals_;
}

void RigidBody::recalcMassMatrix() {
    mass_matrix_ = Matrix(6, 6);
    mass_matrix_(0, 0) = 1.0;
    mass_matrix_(1, 1) = 1.0;
    mass_matrix_(2, 2) = 1.0;
    mass_matrix_ = mass_*mass_matrix_;
    mass_matrix_.Rows(3, 3).Cols(3, 3) = inertia_; // inertia matrix is already multiplied with mass
    mass_matrix_inv_ = inverse(mass_matrix_);
    inertia_inv_ = inverse(inertia_);
  }

VectorView<double> RigidBody::q_trans()
{
  return q_trans_;
}

size_t RigidBody::ConstraintNumberHalf() const
{
  return this->constraint_number_half_;
}

MatrixView<double> RigidBody::q()
{
  return q_;
}

VectorView<double> RigidBody::p_half_trans()
{
  return p_half_trans_;
}

VectorView<double> RigidBody::p_half_skew()
{
  return p_half_skew_;
}

VectorView<double> RigidBody::v_trans()
{
  return v_trans_;
}

VectorView<double> RigidBody::v_skew()
{
  return v_skew_;
}

VectorView<double> RigidBody::p_trans()
{
  return p_trans_;
}

VectorView<double> RigidBody::p_skew()
{
  return p_skew_;
}
MatrixView<double> RigidBody::MassMatrix()
{
  return mass_matrix_;
}
MatrixView<double> RigidBody::MassMatrixInv()
{
  return mass_matrix_inv_;
}
Matrix<double> RigidBody::TensorD()
{
  double tr=0.5*(inertia_(0,0)+inertia_(1,1)+inertia_(2,2));
  Matrix<double> D(3,3);
  for(size_t i=0;i<3;i++) for(size_t j=0;j<3;j++) D(i,j)=(i==j?tr:0.0)-inertia_(i,j);
  return D;
}
MatrixView<double> RigidBody::InertiaInv()
{
  return inertia_inv_;
}
MatrixView<double> RigidBody::Inertia()
{
  return inertia_;
}
double RigidBody::Mass()
{
  return mass_;
}

double& RigidBody::lambda()
{
  return lambda_;
}

double& RigidBody::mu()
{
  return mu_;
}
Vector<double> RigidBody::Vector_q()
{
  Vector q_vec(12);
  q_vec.segment(0, 3) = q_trans_;
  q_vec.segment(3, 9) = ToVector(q_);
  return q_vec;
}

size_t& RigidBody::Index()
{
  return this->body_index_;
}

std::vector<size_t>& RigidBody::Beams()
{
  return this->beams_;
}

void RigidBody::addBeam(size_t i, size_t constraint_number)
{
  if (std::find(this->beams_.begin(), this->beams_.end(), i) == this->beams_.end())
  {
    this->constraint_number_half_ += constraint_number;
    this->beams_.push_back(i);
    this->constraints_.reset(new Matrix(this->constraint_number_half_, 12));
    this->constraints_->setConstant(0);
  }
}

Matrix<double>& RigidBody::Constraints()
{
  return *(this->constraints_);
}

Vector<double>& RigidBody::Force()
{
  return this->force_;
}

// Add a vertex
    void RigidBody::add_vertex(const std::array<double, 3>& v) {
        vertices_.push_back(v);
    }

    // Add a normal
    void RigidBody::add_normal(const std::array<double, 3>& n) {
        normals_.push_back(n);
    }

    // Get vertices
    std::vector<std::array<double, 3>> RigidBody::vertices() {
        return vertices_;
    }

    // Get normals
    std::vector<std::array<double, 3>> RigidBody::normals() {
        return normals_;
    }

#endif