#ifndef EQUATION_RBS_CC
#define EQUATION_RBS_CC

#include "equation_rbs.h"
#include <chrono>
#include "chunked_jacobian.h"

RigidBodySystemEquation::RigidBodySystemEquation(RigidBodySystem& rbs, double step_size): x_(rbs.NumBodies()*dim_per_body() + rbs.NumConstraints()), rbs_(std::ref(rbs)), h_(step_size)
{
  this->dimF_ = rbs.NumBodies()*dim_per_body() + rbs.NumConstraints();
  this->dimX_ = rbs.NumBodies()*dim_per_body() + rbs.NumConstraints();

  this->x_ = rbs.ExpandState();
  rbs.ManageConstraints(x_);
  Vector<double> f(30);
  for (size_t i = 0; i < rbs.NumBodies(); i++)
  {
    std::shared_ptr<RigidBodyEquation> eq = std::make_shared<RigidBodyEquation>(rbs.Bodies(i), rbs, step_size);
    this->func_.addFunction(eq);
  }

  size_t index_count = rbs.BodyDimensions();
  for (size_t i = 0; i < rbs.NumBeams(); i++)
  {
    Beam& bm = rbs.Beams(i);
    bm.StartIndex() = index_count;
    index_count += bm.NumberOfConstraints();
  }
}

size_t RigidBodySystemEquation::DimX() const {
  return this->dimX_;
}
size_t RigidBodySystemEquation::DimF() const
{
  return this->dimF_;
}
size_t RigidBodySystemEquation::BodyDimensions() const
{
  return this->rbs_.get().NumBodies()*dim_per_body();
}

Vector<double>& RigidBodySystemEquation::x()
{
  return this->x_;
}

void RigidBodySystemEquation::Evaluate (VectorView<double> x, VectorView<double> f) const
{
  func_.Evaluate(x, f.segment(0, this->BodyDimensions()));
  size_t constraint_counter = this->BodyDimensions();
  for (size_t i = 0; i < rbs_.get().NumBeams(); i++ )
  {
    Beam& bm = rbs_.get().Beams(i);
    f.segment(constraint_counter, bm.NumberOfConstraints()/2) = this->rbs_.get().Constraint(x, i);
    f.segment(constraint_counter + bm.NumberOfConstraints()/2, bm.NumberOfConstraints()/2) = this->rbs_.get().Velocity_Constraint(x, i);

    constraint_counter += bm.NumberOfConstraints();
  }
}

void RigidBodySystemEquation::EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const
{
  df.setConstant(0);

  const size_t body_rows = this->BodyDimensions();
  func_.EvaluateDeriv(x, df.Block(0, 0, body_rows, dimX_));

  // Beam i's rows read only q (0..11) and the momenta (18..23) of its free bodies.
  RigidBodySystem& rbs = rbs_.get();
  size_t row = body_rows;
  for (size_t i = 0; i < rbs.NumBeams(); i++)
  {
    Beam& bm = rbs.Beams(i);
    const size_t half = bm.NumberOfConstraints()/2;

    ColumnSet cols;
    auto addBody = [&](Connector& c) {
      if (c.Type() == ConnectorType::FIX) return;
      cols.AddRange(dim_per_body()*c.BodyIndex(), 12);
      cols.AddRange(dim_per_body()*c.BodyIndex() + 18, 6);
    };
    addBody(bm.ConnectorA());
    addBody(bm.ConnectorB());

    ChunkedJacobian(x, cols.Finish(), 2*half,
                    [&](auto xd, auto fd) {
                      fd.segment(0, half) = rbs.Constraint(xd, i);
                      fd.segment(half, half) = rbs.Velocity_Constraint(xd, i);
                    },
                    df.Block(row, 0, 2*half, dimX_));
    row += 2*half;
  }
}

void RigidBodySystemEquation::step()
{

  // Start timing
  //auto start_time = std::chrono::high_resolution_clock::now();

  x_ = rbs_.get().ExpandState();

  NewtonSolver((*this), x_, 1e-11, 16);
  rbs_.get().SaveState(x_);

  //auto end_time = std::chrono::high_resolution_clock::now();
  /* std::chrono::duration<double> elapsed = end_time - start_time;
  double elapsed_time = elapsed.count();  // Elapsed time in seconds

  // Sleep if the step was faster than h_
  std::cout << elapsed_time << std::endl;
  if (elapsed_time < 2*h_)
  {
      std::this_thread::sleep_for(std::chrono::duration<double>(h_*2 - elapsed_time));
  } */
  
}

RigidBodySystemEquation assemble(RigidBodySystem& rbs, double step_size)
{
  RigidBodySystemEquation eqrb(rbs, step_size);

  return eqrb;
}

void simulate_rbs(RigidBodySystem& rbs, double step_size, size_t steps_, std::function<void(int,double,VectorView<double>)> callback)

{
  RigidBodySystemEquation eqrb(rbs, step_size);
  
  eqrb.x() = rbs.ExpandState();

  auto start = std::chrono::high_resolution_clock::now();
  for (size_t j = 0; j < steps_; j++) {
    
    NewtonSolver(eqrb, eqrb.x(), 1e-11, 16, callback);
   
    
    rbs.SaveState(eqrb.x());

  
  }
  // Record the end time
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    
}

#endif