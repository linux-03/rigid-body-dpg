#include "equation_rbs.h"

int main()  {
  size_t steps = 400;
  Vector<double> q_(12);
  q_.setConstant(0);
  q_(3) = 1;
  q_(11) = 1;
  q_(7) = 1;
  Vector<double> p_(6);
  p_.setConstant(0);
  //q_(1) = 1;
  q_(2) = 1;

  RigidBodySystem rbs;

  RigidBody rb(q_, p_);
  q_(2) = 2;
  RigidBody rb1(q_, p_);

  rbs.add(rb);
  rbs.add(rb1);

  Connector connect_a({0, 0, 0});
  Connector connect_b({0, 0, 0}, rb, ConnectorType::FREE);
  Connector connector_c({0,0,0}, rb1, ConnectorType::FREE);


  Beam bm(connect_a, connect_b);
  Beam bm1(connect_b, connector_c);

  rbs.add(bm);
  rbs.add(bm1);

  RigidBodySystemEquation eq = assemble(rbs, 0.004);
  for (size_t i = 0; i < steps; i++)
  {
    std::cout << i << std::endl;
    eq.step();
  }
  
  
  /* simulate_rbs(rbs, 0.004, steps, [](int i, double err, VectorView<double> q) {
                      //<<"\t"<< "Translation = " << q(0) << " ," << q(1) << ", "<<", " << q(2) << "} " << std::endl
                      //<<"\t"<< " Rotation: " << q(3) << " ," << q(4) << ", "<<", " << q(5) << "} " << std::endl
                      //<<"\t"<< "           " << q(6) << " ," << q(7) << ", "<<", " << q(8) << "} " << std::endl
                      //<<"\t"<< "           " << q(9) << " ," << q(10) << ", "<<", " << q(11) << "} " << std::endl << std::endl;
                      }
                  ); 
  //}  */

  
  for (size_t j = 0; j < steps/10; j++)
  {
    simulate_rbs(rbs, 0.004, 10, [](int i, double err, VectorView<double> q) {
                      //<<"\t"<< "Translation = " << q(0) << " ," << q(1) << ", "<<", " << q(2) << "} " << std::endl
                      //<<"\t"<< " Rotation: " << q(3) << " ," << q(4) << ", "<<", " << q(5) << "} " << std::endl
                      //<<"\t"<< "           " << q(6) << " ," << q(7) << ", "<<", " << q(8) << "} " << std::endl
                      //<<"\t"<< "           " << q(9) << " ," << q(10) << ", "<<", " << q(11) << "} " << std::endl << std::endl;
                      }
                  ); 
  }
  

  return 0;
}

