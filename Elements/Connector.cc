#ifndef CONNECTOR_CC
#define CONNECTOR_CC

#include "Connector.h"

Connector::Connector(Vector<double> pos, RigidBody& rb, ConnectorType type)
{
  this->type_ = type;
  this->body_index_ = rb.Index();
  this->pos_ = pos;
  this->initial_body_rot_ = rb.q();
  this->initial_body_trans_ = rb.q_trans();
}
Connector::Connector(Vector<double> pos)
{
  this->type_ = ConnectorType::FIX;
  this->body_index_ = 0;
  this->pos_ = pos;
}
Connector::Connector(const Connector& other)
{
  this->type_ = other.type_;
  this->body_index_ = other.body_index_;
  this->pos_ = other.pos_;
  this->initial_body_rot_ = other.initial_body_rot_;
  this->initial_body_trans_ = other.initial_body_trans_;
}
size_t& Connector::BodyIndex()
{
  return this->body_index_;
}
ConnectorType& Connector::Type()
{
  return this->type_;
}
Vector<double> Connector::RefPosition()
{
  return pos_;
}
double Connector::RefPosition(size_t i)
{
  return pos_(i);
}
Vector<double> Connector::InitialPosition()
{
  if (this->type_ == ConnectorType::FIX) return this->pos_;
  return this->initial_body_trans_ + this->initial_body_rot_ * this->pos_;
}
Vector<double> Connector::InitialBodyTranslation()
{
  return this->initial_body_trans_;
}
Matrix<double> Connector::InitialBodyRotation()
{
  return this->initial_body_rot_;
}

#endif