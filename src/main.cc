#include "matrix.h"
#include "vector.h"
#include "rigid_body_system.h"
#include "equation_rbs.h"

int main()  {
  size_t steps = 100;
  Vector<double> q_(12);
  q_.setConstant(0);
  q_(3) = 1;
  q_(11) = 1;
  q_(7) = 1;
  Vector<double> p_(6);
  p_.setConstant(0);
  q_(2) = 2;

  RigidBodySystem rbs;

  RigidBody rb(q_, p_);

  q_(2) = 4;

  RigidBody rb2(q_, p_);

  rbs.add(rb);
  rbs.add(rb2);

  Connector connect_a({0, 0, 0});
  Connector connect_b({0, 0, 0}, rb, ConnectorType::SPHERICAL);
  Connector connect_c({0, 0, 0}, rb, ConnectorType::FREE);
  Connector connect_d({0, 0, 0}, rb2, ConnectorType::FREE);

  Beam bm(connect_a, connect_b);
  Beam bm2(connect_c, connect_d);

  rbs.add(bm);
  rbs.add(bm2);

  //RigidBodySystemEquation eqrb(rbs, 0.1);
  simulate_rbs(rbs, 0.04, 1, [](int i, double err, VectorView<double> q) {
                    std::cout << std::scientific << "Body1 newton-iteration: " << i << " newton-error: " << err << std::endl
                    <<"\t"<< "Translation = " << q(0) << " ," << q(1) << ", "<<", " << q(2) << "} " << std::endl
                    <<"\t"<< " Rotation: " << q(3) << " ," << q(4) << ", "<<", " << q(5) << "} " << std::endl
                    <<"\t"<< "           " << q(6) << " ," << q(7) << ", "<<", " << q(8) << "} " << std::endl
                    <<"\t"<< "           " << q(9) << " ," << q(10) << ", "<<", " << q(11) << "} " << std::endl << std::endl;
                    std::cout << std::scientific << "Body1 newton-iteration: " << i << " newton-error: " << err << std::endl
                    <<"\t"<< "Translation = " << q(0) << " ," << q(1) << ", "<<", " << q(2) << "} " << std::endl
                    <<"\t"<< " Rotation: " << q(3) << " ," << q(4) << ", "<<", " << q(5) << "} " << std::endl
                    <<"\t"<< "           " << q(6) << " ," << q(7) << ", "<<", " << q(8) << "} " << std::endl
                    <<"\t"<< "           " << q(9) << " ," << q(10) << ", "<<", " << q(11) << "} " << std::endl << std::endl;
                      }
                  );

  return 0;
}