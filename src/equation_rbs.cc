#ifndef EQUATION_RBS_CC
#define EQUATION_RBS_CC

#include "equation_rbs.h"
#include <chrono>

RigidBodySystemEquation::RigidBodySystemEquation(RigidBodySystem& rbs, double step_size): x_(rbs.NumBodies()*dim_per_body() + rbs.NumConstraints()), rbs_(std::ref(rbs)), h_(step_size)
{
  this->dimF_ = rbs.NumBodies()*dim_per_body() + rbs.NumConstraints();
  this->dimX_ = rbs.NumBodies()*dim_per_body() + rbs.NumConstraints();

  this->x_ = rbs.ExpandState();
  //std::cout << "Initial x_: " << x_ << " Size: " << x_.Size() << std::endl;
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
  
  Matrix<double> df_block = df.Block(0, 0,  rbs_.get().NumBodies()*dim_per_body(), dimX_);
  func_.EvaluateDeriv(x, df_block);

  df.Block(0, 0,  rbs_.get().NumBodies()*dim_per_body(), dimX_) = df_block;

  Vector<AutoDiff<30, double>> x_diff = x;
  size_t prev_index = 0;
  Vector<AutoDiff<30, double>> f_diff(rbs_.get().NumConstraints());

  for (size_t bd_index = 0; bd_index < rbs_.get().NumBodies(); bd_index++) {

    for (size_t j = 0; j < 30; j++) {
      x_diff(dim_per_body()*prev_index + j).DValue(j) = 0;
      x_diff(dim_per_body()*bd_index + j).DValue(j) = 1;
    }

    prev_index = bd_index;

    size_t constraint_counter = 0;
    for (size_t i = 0; i < rbs_.get().NumBeams(); i++ )
    {
      Beam& bm = rbs_.get().Beams(i);
      f_diff.segment(constraint_counter, bm.NumberOfConstraints()/2) = this->rbs_.get().Constraint(x_diff, i);
      f_diff.segment(constraint_counter + bm.NumberOfConstraints()/2, bm.NumberOfConstraints()/2) = this->rbs_.get().Velocity_Constraint(x_diff, i);

      constraint_counter += bm.NumberOfConstraints();
    }

    for (size_t j = 0; j < rbs_.get().NumConstraints()/2; j++) {
      df.Row(this->BodyDimensions() + 2*j).segment(dim_per_body()*bd_index, dim_per_body()) = f_diff(2*j); 
      df.Row(this->BodyDimensions() + 2*j + 1).segment(dim_per_body()*bd_index, dim_per_body()) = f_diff(2*j + 1); 
    }
    //std::cout << "Computed df for body " << bd_index << std::endl;
  }
  
  
  
  //dNumeric(*(this), x, df);
}

void RigidBodySystemEquation::step()
{

  // Start timing
  //auto start_time = std::chrono::high_resolution_clock::now();

  x_ = rbs_.get().ExpandState();
  //std::cout << "x_: " << x_ << std::endl;

  NewtonSolver((*this), x_, 1e-11, 16);
  //std::cout << x_.segment(3, 9) << std::endl;
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
  // std::cout << "Initial State: " << eqrb.x() << std::endl;

  auto start = std::chrono::high_resolution_clock::now();
  for (size_t j = 0; j < steps_; j++) {
    
    NewtonSolver(eqrb, eqrb.x(), 1e-11, 16, callback);
   
    //std::cout << "step: " << i << std::endl;
    
    //std::cout << eqrb.x().segment(3, 9) << std::endl;
    rbs.SaveState(eqrb.x());
    //x(30)=0;
    //x(31)=0;

    //std::cout << x << std::endl;
    //std::cout << i << ": " << x(1)*9.81 + x.segment(18, 3).squaredNorm() << std::endl << std::endl; 
    //Energy : std::cout << x(1)*9.81 + 0.5*x.segment(18, 3).squaredNorm() << std::endl;
    // Pot Energy : std::cout << x(1)*9.81 << std::endl;
    // std::cout << 0.5*x.segment(18, 3).squaredNorm() << std::endl;
    //std::cout << x(2) << std::endl;
  
  }
  // Record the end time
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    //std::cout << "Time taken: " << duration.count() << " microseconds" << std::endl;
    
}

#endif