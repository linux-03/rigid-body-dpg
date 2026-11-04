#ifndef SPRING_H
#define SPRING_H

#include "Connector.h"
#include <vector.h>

#include <cmath>

/// Linear spring between two connectors, with potential 1/2 * k * (|d| - L)^2,
/// where d is the vector between the connectors and L the rest length.
class Spring{
    size_t index_;
    double stiffness_;
    double rest_length_;
    Connector connector_a_;
    Connector connector_b_;

  public:
    /// A negative rest length means: use the distance when added to the system.
    Spring (Connector connector_a, Connector connector_b, double stiffness = 1.0, double rest_length = -1.0);

    size_t& Index();
    double& Stiffness();
    double& RestLength();
    Connector& ConnectorA();
    Connector& ConnectorB();

    /// Force the spring applies at connector A; connector B receives the opposite force.
    template<typename T>
    Vector<T> ForceOnA(const VectorView<T> q_a, const VectorView<T> q_b) {
      using std::sqrt;
      Vector<T> d = connector_a_.AbsPosition(q_a) - connector_b_.AbsPosition(q_b);
      Vector<T> force(3);
      if (rest_length_ == 0.0) {
        // No division here, so the force stays defined when both ends coincide.
        for (size_t i = 0; i < 3; i++) force(i) = (-stiffness_)*d(i);
        return force;
      }
      T length_squared = d(0)*d(0) + d(1)*d(1) + d(2)*d(2);
      T length = sqrt(length_squared);
      T scale = (stiffness_*rest_length_)/length;
      for (size_t i = 0; i < 3; i++) force(i) = scale*d(i) + (-stiffness_)*d(i);
      return force;
    }

    double Energy(const VectorView<double> q_a, const VectorView<double> q_b) {
      Vector<double> d = connector_a_.AbsPosition(q_a) - connector_b_.AbsPosition(q_b);
      const double stretch = std::sqrt(d(0)*d(0) + d(1)*d(1) + d(2)*d(2)) - rest_length_;
      return 0.5*stiffness_*stretch*stretch;
    }
};

#endif
