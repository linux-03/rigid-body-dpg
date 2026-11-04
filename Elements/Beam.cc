#ifndef CONNECTOR_CC
#define CONNECTOR_CC

#include "Beam.h"
#include "rigid_body.h"

Beam::Beam(Connector connector_a, Connector connector_b): lambda_vec_(2), mu_vec_(2)
{ 
  Vector<double> axis = {0, 0, 0};
  switch(connector_a.Type())
  {
    case ConnectorType::FREE:
      this->number_of_constraints_ += 0;
      break;
    case ConnectorType::SPHERICAL:
      this->number_of_constraints_ += 4;
      axis = connector_a.InitialBodyRotation().transpose() * (connector_b.InitialPosition() - connector_a.InitialBodyTranslation());
      
      this->axis_a_ = axis;
      this->constraint_plane_a_ = normal_plane_basis(axis);
      
      break;
    case ConnectorType::FIX:
      break;
  }
  switch(connector_b.Type())
  {
    case ConnectorType::FREE:
      this->number_of_constraints_ += 0;
      break;
    case ConnectorType::SPHERICAL:

      this->number_of_constraints_ += 4;
      axis = connector_b.InitialBodyRotation().transpose() * (connector_a.InitialPosition() - connector_b.InitialBodyTranslation());
      
      this->axis_b_ = axis;
      this->constraint_plane_b_ = normal_plane_basis(axis);

      break;
    case ConnectorType::FIX:
      break;
  }
  //std::cout << "axis: " << axis << std::endl;
  //std::cout << "Number of constraints in beam: " << this->number_of_constraints_ << std::endl;
  this->lambda_vec_ = Vector<double>(this->number_of_constraints_/2);
  this->mu_vec_ = Vector<double>(this->number_of_constraints_/2);
  this->connector_a_ = connector_a;
  this->connector_b_ = connector_b;
  this->length_ = 0;
  this->index_ = -1;
  //std::cout << "Number of constraints in beam: " << this->number_of_constraints_ << std::endl;
}

Beam::Beam(const Beam& other): lambda_vec_(other.lambda_vec_.Size()), mu_vec_(other.mu_vec_.Size())
{
  this->lambda_vec_ = other.lambda_vec_;
  this->mu_vec_ = other.mu_vec_;
  this->connector_a_ = other.connector_a_;
  this->connector_b_ = other.connector_b_;
  this->number_of_constraints_ = other.number_of_constraints_;
  this->start_index_ = other.start_index_;
  this->constraint_plane_b_ = other.constraint_plane_b_;
  this->constraint_plane_a_ = other.constraint_plane_a_;
  this->length_ = other.length_;
  this->index_ = other.index_;
  this->axis_b_ = other.axis_b_;
  this->axis_a_ = other.axis_a_;
}

size_t& Beam::Index()
{
  return this->index_;
}
double& Beam::Length()
{
  return this->length_;
}
size_t Beam::BodyIndexA()
{
  return this->connector_a_.BodyIndex();
}
size_t Beam::BodyIndexB()
{
  return this->connector_b_.BodyIndex();
}
size_t& Beam::StartIndex()
{
  return this->start_index_;
}
Connector& Beam::ConnectorA()
{
  return this->connector_a_;
}
Connector& Beam::ConnectorB()
{
  return connector_b_;
}
Vector<double>& Beam::LambdaVec()
{
  return this->lambda_vec_;
}
Vector<double>& Beam::MuVec()
{
  return this->mu_vec_;
}
size_t Beam::NumberOfConstraints() {
  return this->number_of_constraints_;
}
#endif