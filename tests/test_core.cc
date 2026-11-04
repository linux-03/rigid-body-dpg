// Regression tests for the rigid body integrator.
//
// Run all:  ctest --test-dir <build> --output-on-failure
// Run one:  <build>/test_core <test_name>
//
// State convention: p_trans = m*v is the world-frame linear momentum and
// p_skew = 1/2*J*omega is the body-frame momentum (see thesis eq. 4.13).

#include "equation_rbs.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void check(bool ok, const std::string& what) {
  if (!ok) throw std::runtime_error(what);
}

void checkNear(double got, double want, double tol, const std::string& what) {
  if (!(std::fabs(got - want) <= tol)) {
    throw std::runtime_error(what + ": got " + std::to_string(got) + ", want " +
                             std::to_string(want) + " +- " + std::to_string(tol));
  }
}

Vector<double> restState(double x, double y, double z) {
  Vector<double> q(12);
  q.setConstant(0);
  q(0) = x; q(1) = y; q(2) = z;
  q(3) = 1; q(7) = 1; q(11) = 1;
  return q;
}

Vector<double> zeroMomentum() {
  Vector<double> p(6);
  p.setConstant(0);
  return p;
}

Matrix<double> diagonalInertia(double a, double b, double c) {
  Matrix<double> J(3, 3);
  J(0, 0) = a; J(1, 1) = b; J(2, 2) = c;
  return J;
}

// ---------------------------------------------------------------------------
// Heavy top benchmark, thesis table 1
// ---------------------------------------------------------------------------
const double kTopMass = 15.0;
const double kTopJ[3] = {0.234375, 0.46875, 0.234375};
const double kTopW0[3] = {0.0, 150.0, -4.61538};
const double kG = 9.81;

struct TopState { double R[9]; double w[3]; };

double crossComponent(const double* a, const double* b, int i) {
  int j = (i + 1) % 3, k = (i + 2) % 3;
  return a[j] * b[k] - a[k] * b[j];
}

// Independent reference: Euler's equations about the fixed point, unrelated to
// the Hamilton-Pontryagin scheme under test.
TopState topRate(const TopState& y) {
  TopState o;
  const double Jo[3] = {kTopJ[0] + kTopMass, kTopJ[1], kTopJ[2] + kTopMass};
  const double xb[3] = {0, 1, 0};
  double gb[3], Jw[3];
  for (int i = 0; i < 3; i++) gb[i] = kTopMass * y.R[6 + i] * (-kG);
  for (int i = 0; i < 3; i++) Jw[i] = Jo[i] * y.w[i];
  for (int i = 0; i < 3; i++) o.w[i] = (crossComponent(xb, gb, i) - crossComponent(y.w, Jw, i)) / Jo[i];
  const double W[9] = {0, -y.w[2], y.w[1], y.w[2], 0, -y.w[0], -y.w[1], y.w[0], 0};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double s = 0;
      for (int k = 0; k < 3; k++) s += y.R[3 * i + k] * W[3 * k + j];
      o.R[3 * i + j] = s;
    }
  return o;
}

TopState addScaled(const TopState& a, const TopState& b, double c) {
  TopState o;
  for (int i = 0; i < 9; i++) o.R[i] = a.R[i] + c * b.R[i];
  for (int i = 0; i < 3; i++) o.w[i] = a.w[i] + c * b.w[i];
  return o;
}

// COM trajectory on a uniform grid of `intervals` cells over [0, 1].
std::vector<std::array<double, 3>> topReference(int intervals) {
  TopState y;
  for (int i = 0; i < 9; i++) y.R[i] = (i % 4 == 0) ? 1.0 : 0.0;
  for (int i = 0; i < 3; i++) y.w[i] = kTopW0[i];
  const int sub = 20;
  const double dt = 1.0 / intervals / sub;
  std::vector<std::array<double, 3>> ref(intervals + 1);
  for (int n = 0; n <= intervals; n++) {
    ref[n] = {y.R[1], y.R[4], y.R[7]};
    if (n == intervals) break;
    for (int k = 0; k < sub; k++) {
      TopState k1 = topRate(y), k2 = topRate(addScaled(y, k1, .5 * dt));
      TopState k3 = topRate(addScaled(y, k2, .5 * dt)), k4 = topRate(addScaled(y, k3, dt));
      TopState t = addScaled(y, k1, dt / 6);
      t = addScaled(t, k2, dt / 3);
      t = addScaled(t, k3, dt / 3);
      y = addScaled(t, k4, dt / 6);
    }
  }
  return ref;
}

struct HeavyTop {
  RigidBodySystem sys;
  Matrix<double> J;
  RigidBody body;
  Connector anchor;
  Connector joint;
  Beam beam;
  Vector<double> offset;

  static RigidBody& placed(RigidBody& b, const Vector<double>& off) {
    b.q_trans()(0) = off(0);
    b.q_trans()(1) = off(1) + 1.0;
    b.q_trans()(2) = off(2);
    for (int i = 0; i < 3; i++) b.p_skew()(i) = 0.5 * kTopJ[i] * kTopW0[i];
    b.p_trans()(0) = kTopMass * 4.61538;
    return b;
  }

  // The body must be registered before its connector reads rb.Index().
  static RigidBody& registered(RigidBodySystem& s, RigidBody& b) { s.add(b); return b; }

  explicit HeavyTop(Vector<double> off = Vector<double>{0, 0, 0})
      : J(diagonalInertia(kTopJ[0], kTopJ[1], kTopJ[2])),
        body(kTopMass, J),
        anchor(off),
        joint(Vector<double>{0, 0, 0}, registered(sys, placed(body, off)), ConnectorType::SPHERICAL),
        beam(anchor, joint),
        offset(off) {
    sys.Gravity().setConstant(0);
    sys.Gravity()(2) = -kG;
    sys.add(beam);
  }

  double energy() {
    double E = kTopMass * kG * body.q_trans()(2);
    for (int i = 0; i < 3; i++) {
      E += 0.5 * body.p_trans()(i) * body.p_trans()(i) / kTopMass;
      const double w = 2.0 * body.p_skew()(i) / kTopJ[i];
      E += 0.5 * kTopJ[i] * w * w;
    }
    return E;
  }
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

// Velocity is p/m: the same physical motion must result for any mass.
void test_mass_scales_kinematics() {
  double x[2], y[2];
  const double masses[2] = {1.0, 15.0};
  for (int i = 0; i < 2; i++) {
    const double m = masses[i];
    Matrix<double> J = diagonalInertia(1, 1, 1);
    RigidBodySystem s;
    RigidBody b(m, J);
    b.q_trans()(2) = 1;
    b.p_trans()(0) = m * 2.0;
    s.add(b);
    Connector a(Vector<double>{0, 0, 0});
    Connector c(Vector<double>{0, 0, 0}, b, ConnectorType::FREE);
    Beam bm(a, c);
    s.add(bm);
    simulate_rbs(s, 0.001, 100);
    x[i] = b.q_trans()(0);
    y[i] = b.q_trans()(1);
  }
  checkNear(x[1], x[0], 1e-9, "x depends on mass");
  checkNear(y[1], y[0], 1e-9, "y depends on mass");
}

// Euler's equations for a torque-free asymmetric body.
void test_inertia_torque_free_precession() {
  const double J[3] = {1, 2, 3};
  double W[3] = {0.1, 0.3, 0.1};
  Matrix<double> Jm = diagonalInertia(J[0], J[1], J[2]);
  RigidBodySystem s;
  s.Gravity().setConstant(0);
  RigidBody b(1.0, Jm);
  for (int i = 0; i < 3; i++) b.p_skew()(i) = 0.5 * J[i] * W[i];
  s.add(b);

  const int N = 20000;
  const double dt = 1.0 / N;
  auto rate = [&](const double* w, double* o) {
    o[0] = (J[1] - J[2]) * w[1] * w[2] / J[0];
    o[1] = (J[2] - J[0]) * w[2] * w[0] / J[1];
    o[2] = (J[0] - J[1]) * w[0] * w[1] / J[2];
  };
  for (int n = 0; n < N; n++) {
    double k1[3], k2[3], k3[3], k4[3], y[3];
    rate(W, k1);
    for (int i = 0; i < 3; i++) y[i] = W[i] + .5 * dt * k1[i];
    rate(y, k2);
    for (int i = 0; i < 3; i++) y[i] = W[i] + .5 * dt * k2[i];
    rate(y, k3);
    for (int i = 0; i < 3; i++) y[i] = W[i] + dt * k3[i];
    rate(y, k4);
    for (int i = 0; i < 3; i++) W[i] += dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
  }

  simulate_rbs(s, 0.001, 1000);
  for (int i = 0; i < 3; i++)
    checkNear(2.0 * b.p_skew()(i), J[i] * W[i], 1e-5, "body angular momentum");
}

// Bodies fall down at the same rate whatever their mass, in a chosen direction.
void test_gravity_direction_and_equivalence() {
  const double masses[2] = {1.0, 15.0};
  for (double m : masses) {
    Matrix<double> J = diagonalInertia(1, 1, 1);
    RigidBodySystem s;
    RigidBody b(m, J);
    s.add(b);
    simulate_rbs(s, 0.001, 100);
    checkNear(b.q_trans()(1), -0.04905, 1e-9, "default gravity fall");

    RigidBodySystem z;
    z.Gravity().setConstant(0);
    z.Gravity()(2) = -9.81;
    RigidBody bz(m, J);
    z.add(bz);
    simulate_rbs(z, 0.001, 100);
    checkNear(bz.q_trans()(1), 0.0, 1e-12, "configured gravity leaks into y");
    checkNear(bz.q_trans()(2), -0.04905, 1e-9, "configured gravity fall");
  }
}

// PotentialGradient holds the force, so it must equal -dV/dq.
void test_potential_consistent_with_force() {
  Matrix<double> J = diagonalInertia(1, 1, 1);
  RigidBodySystem s;
  RigidBody b(15.0, J);
  s.add(b);
  Vector<double> x = s.ExpandState();
  const double eps = 1e-6;
  Vector<double> force = s.PotentialGradient<double>(x, 0);
  for (int i = 0; i < 3; i++) {
    Vector<double> up = x, down = x;
    up(i) += eps;
    down(i) -= eps;
    const double dV = (s.Potential<double>(up) - s.Potential<double>(down)) / (2 * eps);
    checkNear(force(i), -dV, 1e-6, "force is not -dV/dq");
  }
}

// Hidden constraint G*M^-1*p = 0: the COM must move as omega x r on the top.
void test_hidden_constraint_uses_mass() {
  HeavyTop top;
  simulate_rbs(top.sys, 0.0005, 10);
  Matrix<double> R(3, 3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) R(i, j) = top.body.q()(i, j);
  const double wb[3] = {2 * top.body.p_skew()(0) / kTopJ[0],
                        2 * top.body.p_skew()(1) / kTopJ[1],
                        2 * top.body.p_skew()(2) / kTopJ[2]};
  const double xb[3] = {0, 1, 0};
  double vb[3];
  for (int i = 0; i < 3; i++) vb[i] = crossComponent(wb, xb, i);
  for (int i = 0; i < 3; i++) {
    double v = 0;
    for (int j = 0; j < 3; j++) v += R(i, j) * vb[j];
    checkNear(top.body.p_trans()(i) / kTopMass, v, 1e-6, "COM velocity != omega x r");
  }
}

// A fixed connector must not couple the anchored body to body 0.
void test_anchor_on_body_other_than_zero() {
  RigidBodySystem s;
  Vector<double> q0 = restState(1, 0, 0), q1 = restState(0, 0, 2), p = zeroMomentum();
  RigidBody b0(q0, p), b1(q1, p);
  s.add(b0);
  s.add(b1);
  Connector a0(Vector<double>{0, 0, 0}), a1(Vector<double>{0, 0, 0});
  Connector c0(Vector<double>{0, 0, 0}, b0, ConnectorType::FREE);
  Connector c1(Vector<double>{0, 0, 0}, b1, ConnectorType::FREE);
  Beam beam0(a0, c0), beam1(a1, c1);
  s.add(beam0);
  s.add(beam1);
  b0.p_trans()(2) = 1;
  b1.p_trans()(0) = 1;

  Vector<double> x = s.ExpandState();
  checkNear(s.Velocity_Constraint<double>(x, 1)(0), 0.0, 1e-12, "velocity constraint at start");

  simulate_rbs(s, 0.001, 1);
  // h*v^2/(2r) = 0.0005 is the true centripetal change; the old bug gave about 1.
  checkNear(b1.p_trans()(2), 0.0, 1e-3, "spurious radial momentum on body 1");
}

// The kinematic equation must not be singular at exactly representable rotations.
void test_exact_rotations_are_regular() {
  for (int mode = 0; mode < 4; mode++) {
    RigidBodySystem s;
    Vector<double> q = restState(0, 0, 1), p = zeroMomentum();
    RigidBody b(q, p);
    MatrixView<double> Q = b.q();
    Q.setConstant(0);
    if (mode == 0) { Q(0, 1) = -1; Q(1, 0) = 1; Q(2, 2) = 1; }   // Rz(90)
    if (mode == 1) { Q(0, 0) = 1; Q(1, 1) = -1; Q(2, 2) = -1; }  // Rx(180)
    if (mode == 2) { Q(0, 0) = -1; Q(1, 1) = -1; Q(2, 2) = 1; }  // Rz(180)
    if (mode == 3) { Q(0, 2) = 1; Q(1, 1) = 1; Q(2, 0) = -1; }   // Ry(90)
    b.p_skew()(0) = 0.5; b.p_skew()(1) = 1; b.p_skew()(2) = 1.5;
    s.add(b);
    Connector a(Vector<double>{0, 0, 0});
    Connector c(Vector<double>{0, 0, 0}, b, ConnectorType::FREE);
    Beam bm(a, c);
    s.add(bm);
    simulate_rbs(s, 0.01, 5);
  }
}

// The dynamics must not depend on the world orientation of a free body.
void test_left_invariance() {
  double reference[3] = {0, 0, 0};
  for (int mode = 0; mode < 2; mode++) {
    Matrix<double> J = diagonalInertia(1, 2, 3);
    RigidBodySystem s;
    s.Gravity().setConstant(0);
    RigidBody b(1.0, J);
    if (mode == 1) {
      MatrixView<double> Q = b.q();
      Q.setConstant(0);
      Q(0, 2) = 1; Q(1, 1) = 1; Q(2, 0) = -1;
    }
    b.p_skew()(0) = 0.05; b.p_skew()(1) = 0.3; b.p_skew()(2) = 0.15;
    s.add(b);
    simulate_rbs(s, 0.01, 100);
    for (int i = 0; i < 3; i++) {
      if (mode == 0) reference[i] = b.p_skew()(i);
      else checkNear(b.p_skew()(i), reference[i], 1e-9, "result depends on world orientation");
    }
  }
}

void test_body_without_beams_runs() {
  RigidBodySystem s;
  Vector<double> q = restState(0, 0, 0), p = zeroMomentum();
  RigidBody b(q, p);
  s.add(b);
  simulate_rbs(s, 0.01, 1);
  check(std::isfinite(b.q_trans()(1)), "state is not finite");
}

void test_copy_keeps_index_and_constraints() {
  RigidBodySystem s;
  Vector<double> q = restState(0, 0, 1), p = zeroMomentum();
  RigidBody b(q, p);
  s.add(b);
  Connector a(Vector<double>{0, 0, 0});
  Connector c(Vector<double>{0, 0, 0}, b, ConnectorType::FREE);
  Beam bm(a, c);
  s.add(bm);
  std::vector<RigidBody> copies = s.Bodies();
  check(copies[0].Index() == b.Index(), "copy lost its body index");
  check(copies[0].ConstraintNumberHalf() == b.ConstraintNumberHalf(), "copy lost its constraint count");

  RigidBody lone(q, p);
  RigidBody loneCopy(lone);  // no beams: used to dereference a null pointer
  check(loneCopy.Beams().empty(), "copy of a beamless body has beams");
}

// A fixed anchor away from the origin is the same problem, shifted.
void test_fixed_anchor_away_from_origin() {
  Connector anchor(Vector<double>{1, 0, 0});
  checkNear(anchor.InitialPosition()(0), 1.0, 1e-15, "InitialPosition ignores fixed position");

  HeavyTop origin;
  HeavyTop shifted(Vector<double>{1, 0, 0});
  simulate_rbs(origin.sys, 0.0005, 50);
  simulate_rbs(shifted.sys, 0.0005, 50);
  checkNear(shifted.body.q_trans()(0) - 1.0, origin.body.q_trans()(0), 1e-8, "shifted top x");
  checkNear(shifted.body.q_trans()(1), origin.body.q_trans()(1), 1e-8, "shifted top y");
  checkNear(shifted.body.q_trans()(2), origin.body.q_trans()(2), 1e-8, "shifted top z");
}

// Thesis theorem 4.3: second order, against a reference from a different method.
void test_heavy_top_second_order() {
  const int finest = 4000;
  const auto ref = topReference(finest);
  checkNear(ref[finest][2], -0.748491, 1e-3, "reference solution changed");

  double errors[3];
  const int factors[3] = {4, 2, 1};
  for (int level = 0; level < 3; level++) {
    const int fac = factors[level];
    const double h = 1.0 / finest * fac;
    const int steps = finest / fac;
    HeavyTop top;
    const double E0 = top.energy();
    RigidBodySystemEquation eq(top.sys, h);
    eq.x() = top.sys.ExpandState();
    double err = 0, dE = 0;
    for (int k = 0; k < steps; k++) {
      NewtonSolver(eq, eq.x(), 1e-11, 16);
      top.sys.SaveState(eq.x());
      const auto& r = ref[(k + 1) * fac];
      double e = 0;
      for (int i = 0; i < 3; i++) e += std::pow(top.body.q_trans()(i) - r[i], 2);
      err = std::max(err, std::sqrt(e));
      dE = std::max(dE, std::fabs(top.energy() - E0));
    }
    errors[level] = err;
    check(dE / E0 < 1e-4, "energy drift too large");
    checkNear(std::sqrt(std::pow(top.body.q_trans()(0), 2) + std::pow(top.body.q_trans()(1), 2) +
                        std::pow(top.body.q_trans()(2), 2)), 1.0, 1e-10, "rod length");
  }
  for (int level = 0; level < 2; level++)
    checkNear(std::log2(errors[level] / errors[level + 1]), 2.0, 0.2, "convergence rate");
  check(errors[2] < 1e-4, "heavy top error too large");
}

// Energy of a symplectic method stays bounded over long runs.
void test_pendulum_energy_bounded() {
  RigidBodySystem s;
  Vector<double> q = restState(0, 0, 1), p = zeroMomentum();
  RigidBody b(q, p);
  s.add(b);
  Connector a(Vector<double>{0, 0, 0});
  Connector c(Vector<double>{0, 0, 0}, b, ConnectorType::FREE);
  Beam bm(a, c);
  s.add(bm);
  auto energy = [&]() {
    double E = 9.81 * b.q_trans()(1);
    for (int i = 0; i < 3; i++) E += 0.5 * b.p_trans()(i) * b.p_trans()(i);
    return E;
  };
  const double E0 = energy();
  double early = 0, late = 0;
  const int steps = 20000;
  for (int n = 0; n < steps; n++) {
    simulate_rbs(s, 0.001, 1);
    double d = std::fabs(energy() - E0);
    (n < steps / 2 ? early : late) = std::max(n < steps / 2 ? early : late, d);
  }
  check(early < 1e-3, "energy error too large");
  check(late < 1.1 * early, "energy drifts over time");
}

// With no external forces, energy and angular momentum are constant. Their drift
// measures how exactly each step is solved, so this guards the solver's accuracy.
void test_free_body_invariants_are_tight() {
  const double J[3] = {0.0012, 0.0021, 0.0033};
  const double w[3] = {0.02, 5, 0.02};
  Matrix<double> Jm = diagonalInertia(J[0], J[1], J[2]);
  RigidBodySystem s;
  s.Gravity().setConstant(0);
  RigidBody b(1.0, Jm);
  for (int i = 0; i < 3; i++) b.p_skew()(i) = 0.5 * J[i] * w[i];
  s.add(b);
  auto energy = [&]() {
    double E = 0;
    for (int i = 0; i < 3; i++) {
      const double wi = 2.0 * b.p_skew()(i) / J[i];
      E += 0.5 * J[i] * wi * wi;
    }
    return E;
  };
  auto momentum = [&](int i) {
    double L = 0;
    for (int j = 0; j < 3; j++) L += b.q()(i, j) * 2.0 * b.p_skew()(j);
    return L;
  };
  const double E0 = energy();
  const double L0[3] = {momentum(0), momentum(1), momentum(2)};
  const double L0norm = std::sqrt(L0[0] * L0[0] + L0[1] * L0[1] + L0[2] * L0[2]);

  simulate_rbs(s, 0.004, 2500);

  checkNear(energy() / E0, 1.0, 1e-12, "relative energy");
  for (int i = 0; i < 3; i++)
    checkNear((momentum(i) - L0[i]) / L0norm, 0.0, 1e-12, "relative angular momentum");
}

// Zero rest length gives x'' = -(k/m) x, which has an exact solution.
void test_spring_matches_exact_oscillator() {
  const double m = 2.0, k = 50.0, x0 = 0.3, T = 1.0, omega = std::sqrt(k / m);
  double previous = 0;
  for (int steps : {200, 400}) {
    Matrix<double> J = diagonalInertia(1, 1, 1);
    RigidBodySystem s;
    s.Gravity().setConstant(0);
    RigidBody b(m, J);
    b.q_trans()(0) = x0;
    s.add(b);
    Connector anchor(Vector<double>{0, 0, 0});
    Connector end(Vector<double>{0, 0, 0}, b, ConnectorType::FREE);
    Spring spring(anchor, end, k, 0.0);
    s.add(spring);
    simulate_rbs(s, T / steps, steps);
    const double error = std::fabs(b.q_trans()(0) - x0 * std::cos(omega * T));
    check(error < 2e-3, "oscillator error too large: " + std::to_string(error));
    if (previous > 0) checkNear(std::log2(previous / error), 2.0, 0.1, "oscillator convergence rate");
    previous = error;
  }
}

// The default rest length is the distance when added, so the spring starts relaxed.
void test_spring_default_rest_length() {
  Matrix<double> J = diagonalInertia(1, 1, 1);
  RigidBodySystem s;
  s.Gravity().setConstant(0);
  RigidBody b(1.0, J);
  b.q_trans()(0) = 2;
  s.add(b);
  Connector anchor(Vector<double>{0, 0, 0});
  Connector end(Vector<double>{0, 0, 0}, b, ConnectorType::FREE);
  Spring spring(anchor, end, 10.0);
  s.add(spring);
  checkNear(spring.RestLength(), 2.0, 1e-14, "rest length");
  simulate_rbs(s, 0.01, 10);
  checkNear(b.q_trans()(0), 2.0, 1e-12, "relaxed spring moved the body");
}

// The stored force must equal minus the gradient, including the rotational part.
void test_spring_force_is_minus_gradient() {
  Matrix<double> J = diagonalInertia(1, 2, 3);
  RigidBodySystem s;
  RigidBody a(2.0, J), c(3.0, J);
  a.q_trans()(0) = 1;
  c.q_trans()(1) = 2;
  s.add(a);
  s.add(c);
  Connector end_a(Vector<double>{0.1, 0.2, -0.3}, a, ConnectorType::FREE);
  Connector end_c(Vector<double>{-0.2, 0.1, 0.4}, c, ConnectorType::FREE);
  Spring spring(end_a, end_c, 7.0, 0.5);
  s.add(spring);
  Vector<double> x = s.ExpandState();
  const double eps = 1e-6;
  for (size_t body = 0; body < 2; body++) {
    Vector<double> force = s.PotentialGradient<double>(x, body);
    for (size_t i = 0; i < 12; i++) {
      Vector<double> up = x, down = x;
      up(body * dim_per_body() + i) += eps;
      down(body * dim_per_body() + i) -= eps;
      const double dV = (s.Potential<double>(up) - s.Potential<double>(down)) / (2 * eps);
      checkNear(force(i), -dV, 1e-6, "force is not -dV/dq");
    }
  }
}

// Internal forces leave total linear and angular momentum unchanged.
void test_spring_pair_conserves_momentum() {
  Matrix<double> J = diagonalInertia(1, 2, 3);
  RigidBodySystem s;
  s.Gravity().setConstant(0);
  RigidBody a(2.0, J), c(3.0, J);
  a.q_trans()(0) = 1;
  c.q_trans()(1) = 2;
  a.p_trans()(2) = 1;
  c.p_skew()(0) = 0.4;
  s.add(a);
  s.add(c);
  Connector end_a(Vector<double>{0.1, 0.2, -0.3}, a, ConnectorType::FREE);
  Connector end_c(Vector<double>{-0.2, 0.1, 0.4}, c, ConnectorType::FREE);
  Spring spring(end_a, end_c, 7.0, 0.5);
  s.add(spring);

  auto totals = [&](double* out) {
    RigidBody* bodies[2] = {&a, &c};
    for (int i = 0; i < 6; i++) out[i] = 0;
    for (RigidBody* body : bodies) {
      double q[3], p[3];
      for (int i = 0; i < 3; i++) { q[i] = body->q_trans()(i); p[i] = body->p_trans()(i); }
      for (int i = 0; i < 3; i++) {
        out[i] += p[i];
        out[3 + i] += crossComponent(q, p, i);
        for (int j = 0; j < 3; j++) out[3 + i] += body->q()(i, j) * 2.0 * body->p_skew()(j);
      }
    }
  };
  auto energy = [&]() {
    Vector<double> x = s.ExpandState();
    double E = s.Potential<double>(x);
    RigidBody* bodies[2] = {&a, &c};
    const double masses[2] = {2.0, 3.0}, Jd[3] = {1, 2, 3};
    for (int b = 0; b < 2; b++)
      for (int i = 0; i < 3; i++) {
        E += 0.5 * bodies[b]->p_trans()(i) * bodies[b]->p_trans()(i) / masses[b];
        E += 2.0 * bodies[b]->p_skew()(i) * bodies[b]->p_skew()(i) / Jd[i];
      }
    return E;
  };
  double start[6], end[6];
  totals(start);
  const double E0 = energy();
  double band = 0;
  for (int n = 0; n < 1000; n++) {
    simulate_rbs(s, 0.002, 1);
    band = std::max(band, std::fabs(energy() - E0));
  }
  totals(end);
  for (int i = 0; i < 6; i++) checkNear(end[i], start[i], 1e-9, "total momentum component " + std::to_string(i));
  check(band / E0 < 1e-3, "energy error too large: " + std::to_string(band / E0));
}

void test_spring_jacobian_matches_finite_differences() {
  Matrix<double> J = diagonalInertia(1, 2, 3);
  RigidBodySystem s;
  RigidBody a(2.0, J), c(3.0, J);
  a.q_trans()(0) = 1;
  c.q_trans()(1) = 2;
  s.add(a);
  s.add(c);
  Connector end_a(Vector<double>{0.1, 0.2, -0.3}, a, ConnectorType::FREE);
  Connector end_c(Vector<double>{-0.2, 0.1, 0.4}, c, ConnectorType::FREE);
  Spring spring(end_a, end_c, 7.0, 0.5);
  s.add(spring);
  RigidBodySystemEquation eq(s, 0.001);
  Vector<double> x = s.ExpandState();
  for (size_t i = 0; i < x.Size(); i++) x(i) += 0.01 * std::sin(i + 1.0);
  const size_t n = eq.DimX();
  Matrix<double> exact(n, n), numeric(n, n);
  eq.EvaluateDeriv(x, exact);
  dNumeric(eq, x, numeric);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      checkNear(exact(i, j), numeric(i, j), 1e-5,
                "Jacobian entry (" + std::to_string(i) + "," + std::to_string(j) + ")");
}

// ---------------------------------------------------------------------------
// Symplecticity: a symplectic flow keeps the canonical two-form
//   Omega(u, v) = du_q . dv_p - du_p . dv_q
// constant. We embed a state as (x, Q, p, Q*hat(p_skew)) in R^24, take tangent
// vectors by finite differences, and compare Omega before and after the flow.
// ---------------------------------------------------------------------------
using Embedding = std::array<double, 24>;
using Flow = std::function<Embedding(const std::vector<double>&, int)>;

void rotationFromVector(const double* v, double* R) {
  const double angle = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  const double a = angle > 1e-12 ? std::sin(angle) / angle : 1.0;
  const double b = angle > 1e-6 ? (1.0 - std::cos(angle)) / (angle * angle) : 0.5;
  const double W[9] = {0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double WW = 0;
      for (int k = 0; k < 3; k++) WW += W[3 * i + k] * W[3 * k + j];
      R[3 * i + j] = (i == j ? 1.0 : 0.0) + a * W[3 * i + j] + b * WW;
    }
}

Embedding embed(RigidBody& b) {
  Embedding e;
  const double w[3] = {b.p_skew()(0), b.p_skew()(1), b.p_skew()(2)};
  const double W[9] = {0, -w[2], w[1], w[2], 0, -w[0], -w[1], w[0], 0};
  for (int i = 0; i < 3; i++) {
    e[i] = b.q_trans()(i);
    e[12 + i] = b.p_trans()(i);
    for (int j = 0; j < 3; j++) {
      e[3 + 3 * i + j] = b.q()(i, j);
      double P = 0;
      for (int k = 0; k < 3; k++) P += b.q()(i, k) * W[3 * k + j];
      e[15 + 3 * i + j] = P;
    }
  }
  return e;
}

std::vector<double> twoForm(const Flow& flow, const std::vector<double>& z, int steps) {
  const double eps = 1e-5;
  const size_t n = z.size();
  std::vector<Embedding> tangents(n);
  for (size_t i = 0; i < n; i++) {
    std::vector<double> up = z, down = z;
    up[i] += eps;
    down[i] -= eps;
    const Embedding a = flow(up, steps), b = flow(down, steps);
    for (int k = 0; k < 24; k++) tangents[i][k] = (a[k] - b[k]) / (2 * eps);
  }
  std::vector<double> omega(n * n);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++) {
      double sum = 0;
      for (int k = 0; k < 12; k++)
        sum += tangents[i][k] * tangents[j][12 + k] - tangents[i][12 + k] * tangents[j][k];
      omega[i * n + j] = sum;
    }
  return omega;
}

double symplecticDefect(const Flow& flow, const std::vector<double>& z, int steps) {
  const auto before = twoForm(flow, z, 0), after = twoForm(flow, z, steps);
  double largest = 0, defect = 0;
  for (size_t i = 0; i < before.size(); i++) {
    largest = std::max(largest, std::fabs(before[i]));
    defect = std::max(defect, std::fabs(after[i] - before[i]));
  }
  return defect / largest;
}

// A free body on a nonlinear spring attached off-centre, in 12 unconstrained coordinates.
Embedding springBodyFlow(const std::vector<double>& z, int steps, double damping) {
  const double Jd[3] = {1, 2, 3};
  Matrix<double> J = diagonalInertia(Jd[0], Jd[1], Jd[2]);
  RigidBodySystem s;
  RigidBody b(2.0, J);
  double R[9];
  rotationFromVector(&z[3], R);
  for (int i = 0; i < 3; i++) {
    b.q_trans()(i) = z[i];
    b.p_trans()(i) = z[6 + i];
    b.p_skew()(i) = z[9 + i];
    for (int j = 0; j < 3; j++) b.q()(i, j) = R[3 * i + j];
  }
  s.add(b);
  Connector anchor(Vector<double>{0, 0, 0});
  Connector end(Vector<double>{0.3, 0.1, 0}, b, ConnectorType::FREE);
  Spring spring(anchor, end, 20.0, 0.5);
  s.add(spring);
  for (int n = 0; n < steps; n++) {
    simulate_rbs(s, 0.01, 1);
    for (int i = 0; i < 3; i++) {
      b.p_trans()(i) = damping * b.p_trans()(i);
      b.p_skew()(i) = damping * b.p_skew()(i);
    }
  }
  return embed(b);
}

const std::vector<double> kSpringBodyState = {1, 0.2, -0.3, 0.2, -0.1, 0.3, 0.5, -0.4, 0.3, 0.2, 0.5, -0.3};

void test_symplectic_body_on_spring() {
  Flow flow = [](const std::vector<double>& z, int steps) { return springBodyFlow(z, steps, 1.0); };
  const double defect = symplecticDefect(flow, kSpringBodyState, 50);
  std::printf("symplectic defect, body on spring: %.3e\n", defect);
  check(defect < 1e-7, "two-form changed: " + std::to_string(defect));
}

// Guards the test itself: a slightly damped flow is not symplectic and must be caught.
void test_symplectic_check_detects_damping() {
  Flow flow = [](const std::vector<double>& z, int steps) { return springBodyFlow(z, steps, 0.9999); };
  const double defect = symplecticDefect(flow, kSpringBodyState, 50);
  std::printf("symplectic defect, damped flow: %.3e\n", defect);
  check(defect > 1e-3, "damping was not detected: " + std::to_string(defect));
}

// The constrained heavy top, in 6 coordinates on its manifold: rotation and angular velocity.
void test_symplectic_heavy_top() {
  Flow flow = [](const std::vector<double>& z, int steps) {
    Matrix<double> J = diagonalInertia(kTopJ[0], kTopJ[1], kTopJ[2]);
    RigidBodySystem s;
    s.Gravity().setConstant(0);
    s.Gravity()(2) = -kG;
    RigidBody b(kTopMass, J);
    double R[9];
    rotationFromVector(&z[0], R);
    const double xb[3] = {0, 1, 0}, w[3] = {z[3], z[4], z[5]};
    double v[3];
    for (int i = 0; i < 3; i++) v[i] = crossComponent(w, xb, i);
    for (int i = 0; i < 3; i++) {
      double q = 0, p = 0;
      for (int j = 0; j < 3; j++) {
        b.q()(i, j) = R[3 * i + j];
        q += R[3 * i + j] * xb[j];
        p += R[3 * i + j] * v[j];
      }
      b.q_trans()(i) = q;
      b.p_trans()(i) = kTopMass * p;
      b.p_skew()(i) = 0.5 * kTopJ[i] * w[i];
    }
    s.add(b);
    Connector anchor(Vector<double>{0, 0, 0});
    Connector joint(Vector<double>{0, 0, 0}, b, ConnectorType::SPHERICAL);
    Beam beam(anchor, joint);
    s.add(beam);
    if (steps > 0) simulate_rbs(s, 0.002, steps);
    return embed(b);
  };
  const double defect = symplecticDefect(flow, {0.1, -0.2, 0.3, 1, 20, -2}, 50);
  std::printf("symplectic defect, heavy top: %.3e\n", defect);
  check(defect < 1e-6, "two-form changed: " + std::to_string(defect));
}

// Automatic differentiation must agree with central differences.
void test_jacobian_matches_finite_differences() {
  HeavyTop top;
  RigidBodySystemEquation eq(top.sys, 0.001);
  Vector<double> x = top.sys.ExpandState();
  for (size_t i = 0; i < x.Size(); i++) x(i) += 0.01 * std::sin(i + 1.0);
  const size_t n = eq.DimX();
  Matrix<double> exact(n, n), numeric(n, n);
  eq.EvaluateDeriv(x, exact);
  dNumeric(eq, x, numeric);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      checkNear(exact(i, j), numeric(i, j), 1e-5,
                "Jacobian entry (" + std::to_string(i) + "," + std::to_string(j) + ")");
}

const std::map<std::string, std::function<void()>>& allTests() {
  static const std::map<std::string, std::function<void()>> tests = {
      {"mass_scales_kinematics", test_mass_scales_kinematics},
      {"inertia_torque_free_precession", test_inertia_torque_free_precession},
      {"gravity_direction_and_equivalence", test_gravity_direction_and_equivalence},
      {"potential_consistent_with_force", test_potential_consistent_with_force},
      {"hidden_constraint_uses_mass", test_hidden_constraint_uses_mass},
      {"anchor_on_body_other_than_zero", test_anchor_on_body_other_than_zero},
      {"exact_rotations_are_regular", test_exact_rotations_are_regular},
      {"left_invariance", test_left_invariance},
      {"body_without_beams_runs", test_body_without_beams_runs},
      {"copy_keeps_index_and_constraints", test_copy_keeps_index_and_constraints},
      {"fixed_anchor_away_from_origin", test_fixed_anchor_away_from_origin},
      {"heavy_top_second_order", test_heavy_top_second_order},
      {"pendulum_energy_bounded", test_pendulum_energy_bounded},
      {"free_body_invariants_are_tight", test_free_body_invariants_are_tight},
      {"jacobian_matches_finite_differences", test_jacobian_matches_finite_differences},
      {"spring_matches_exact_oscillator", test_spring_matches_exact_oscillator},
      {"spring_default_rest_length", test_spring_default_rest_length},
      {"spring_force_is_minus_gradient", test_spring_force_is_minus_gradient},
      {"spring_pair_conserves_momentum", test_spring_pair_conserves_momentum},
      {"spring_jacobian_matches_finite_differences", test_spring_jacobian_matches_finite_differences},
      {"symplectic_body_on_spring", test_symplectic_body_on_spring},
      {"symplectic_check_detects_damping", test_symplectic_check_detects_damping},
      {"symplectic_heavy_top", test_symplectic_heavy_top},
  };
  return tests;
}

}  // namespace

int main(int argc, char** argv) {
  const auto& tests = allTests();
  int failed = 0;
  auto run = [&](const std::string& name, const std::function<void()>& fn) {
    try {
      fn();
      std::printf("PASS %s\n", name.c_str());
    } catch (const std::exception& e) {
      std::printf("FAIL %s: %s\n", name.c_str(), e.what());
      failed++;
    }
    std::fflush(stdout);
  };
  if (argc > 1) {
    auto it = tests.find(argv[1]);
    if (it == tests.end()) {
      std::printf("unknown test %s\n", argv[1]);
      return 2;
    }
    run(it->first, it->second);
  } else {
    for (const auto& t : tests) run(t.first, t.second);
  }
  return failed ? 1 : 0;
}
