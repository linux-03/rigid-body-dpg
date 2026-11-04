"""Build a mechanical system, simulate it, and keep the trajectory."""
import time

import numpy as np

from .rigid_body_FEM import (
    Beam,
    Connector,
    ConnectorType,
    RBS_FEM,
    RigidBody_FEM,
    bla,
    simulate,
)

PALETTE = ("#d1495b", "#00798c", "#edae49", "#30638e", "#66a182")


class SimulationError(RuntimeError):
    """The nonlinear solver failed during a time step."""


def _vector(values):
    out = bla.Vector(3)
    for i in range(3):
        out[i] = float(values[i])
    return out


def _matrix(values):
    out = bla.Matrix(3, 3)
    for i in range(3):
        for j in range(3):
            out[i, j] = float(values[i, j])
    return out


def _triangles(solid):
    """Triangle vertices and normals of a netgen solid, as (n, 3) arrays."""
    data = solid._webgui_data()["Bezier_trig_points"]
    vertices = np.stack([np.asarray(data[k], dtype=np.float32).reshape(-1, 4)[:, :3] for k in range(3)], axis=1)
    normals = np.stack([np.asarray(data[k], dtype=np.float32).reshape(-1, 3) for k in range(3, 6)], axis=1)
    return vertices.reshape(-1, 3), normals.reshape(-1, 3)


class Body:
    """One rigid body: its physical data, its shape, and its solver object."""

    def __init__(self, name, mass, inertia, vertices, normals, color, at, velocity, spin):
        self.name = name
        self.mass = float(mass)
        self.inertia = np.asarray(inertia, dtype=float)
        self.inertia_inv = np.linalg.inv(self.inertia)
        self.vertices = vertices
        self.normals = normals
        self.color = color
        self.solver = RigidBody_FEM(self.mass, _matrix(self.inertia))
        self.solver.q_trans = [float(v) for v in at]
        self.solver.p_trans = [float(v) * self.mass for v in velocity]
        # The solver's angular variable is the body-frame momentum J*omega/2.
        self.solver.p_skew = [float(v) for v in 0.5 * self.inertia @ np.asarray(spin, dtype=float)]

    @property
    def position(self):
        return np.array(self.solver.q_trans)

    @property
    def rotation(self):
        m = self.solver.asTuple()
        return np.array([[m[0], m[4], m[8]], [m[1], m[5], m[9]], [m[2], m[6], m[10]]])

    @property
    def angular_velocity(self):
        """Angular velocity in the body frame."""
        return 2.0 * self.inertia_inv @ np.array(self.solver.p_skew)

    @property
    def velocity(self):
        return np.array(self.solver.p_trans) / self.mass

    def kinetic_energy(self):
        w = self.angular_velocity
        v = self.velocity
        return 0.5 * self.mass * v @ v + 0.5 * w @ self.inertia @ w

    def angular_momentum(self):
        """Spin angular momentum about the centre of mass, in the world frame."""
        return self.rotation @ (self.inertia @ self.angular_velocity)


class World:
    """A collection of rigid bodies and joints under uniform gravity."""

    def __init__(self, gravity=(0.0, 0.0, -9.81)):
        self.gravity = np.asarray(gravity, dtype=float)
        self.system = RBS_FEM()
        self.system.gravity = [float(g) for g in self.gravity]
        self.bodies = []
        self.links = []  # (anchor, offset_a, body_b, offset_b); anchor is a Body or a point
        self._keep_alive = []  # the C++ system stores references to these objects

    def add(self, shape, at=None, velocity=(0, 0, 0), spin=(0, 0, 0), density=1.0,
            mass=None, inertia=None, color=None, name=None):
        """Add a body made from a netgen solid.

        The mass and inertia come from the solid and `density`, about its centre
        of mass, unless given explicitly. `at` places the centre of mass and
        defaults to where the solid was modelled. `spin` is the initial angular
        velocity about the body axes.
        """
        center = np.array([shape.center[i] for i in range(3)])
        solid = shape.Move(tuple(-center))
        computed = np.array([[solid.inertia[i, j] for j in range(3)] for i in range(3)])
        vertices, normals = _triangles(solid)
        body = Body(
            name=name or f"body{len(self.bodies)}",
            mass=density * solid.mass if mass is None else mass,
            inertia=density * computed if inertia is None else inertia,
            vertices=vertices,
            normals=normals,
            color=color or PALETTE[len(self.bodies) % len(PALETTE)],
            at=center if at is None else at,
            velocity=velocity,
            spin=spin,
        )
        self.system.add(body.solver)  # must precede any connector on this body
        self.bodies.append(body)
        return body

    def _add_beam(self, connector_a, connector_b):
        beam = Beam(connector_a, connector_b)
        self.system.add(beam)
        self._keep_alive += [connector_a, connector_b, beam]

    @staticmethod
    def _require_distance(p, q, what):
        if np.linalg.norm(np.asarray(p, dtype=float) - np.asarray(q, dtype=float)) < 1e-9:
            raise ValueError(f"{what} needs two distinct points")

    def pin(self, body, point):
        """Ball joint between `body` and a fixed `point`, which stays rigidly attached."""
        point = np.asarray(point, dtype=float)
        self._require_distance(point, body.position, "pin")
        self._add_beam(Connector(_vector(point)),
                       Connector(_vector((0, 0, 0)), body.solver, ConnectorType.SPHERICAL))
        self.links.append((point, None, body, np.zeros(3)))

    def ball(self, a, b, point):
        """Ball joint between bodies `a` and `b` at the world `point`."""
        point = np.asarray(point, dtype=float)
        self._require_distance(point, b.position, "ball")
        offset = a.rotation.T @ (point - a.position)
        self._add_beam(Connector(_vector(offset), a.solver, ConnectorType.FREE),
                       Connector(_vector((0, 0, 0)), b.solver, ConnectorType.SPHERICAL))
        self.links.append((a, offset, b, np.zeros(3)))

    def rod(self, a, b, at_a=(0, 0, 0), at_b=(0, 0, 0)):
        """Massless rod of fixed length; both ends pivot freely. Offsets are body-frame."""
        at_a, at_b = np.asarray(at_a, dtype=float), np.asarray(at_b, dtype=float)
        self._require_distance(a.position + a.rotation @ at_a, b.position + b.rotation @ at_b, "rod")
        self._add_beam(Connector(_vector(at_a), a.solver, ConnectorType.FREE),
                       Connector(_vector(at_b), b.solver, ConnectorType.FREE))
        self.links.append((a, at_a, b, at_b))

    def potential_energy(self):
        return -sum(b.mass * self.gravity @ b.position for b in self.bodies)

    def energy(self):
        return sum(b.kinetic_energy() for b in self.bodies) + self.potential_energy()

    def _link_ends(self):
        ends = []
        for a, offset_a, b, offset_b in self.links:
            start = a if offset_a is None else a.position + a.rotation @ offset_a
            ends.append((start, b.position + b.rotation @ offset_b))
        return ends

    def run(self, duration, h=0.002, fps=30):
        """Advance the system by `duration` seconds with time step `h`.

        Returns a Recording sampled at about `fps` frames per simulated second.
        """
        per_frame = max(1, round(1.0 / (fps * h)))
        total = max(1, round(duration / h))
        rec = Recording(self)
        rec._capture(0.0)
        started = time.perf_counter()
        for done in range(0, total, per_frame):
            count = min(per_frame, total - done)
            try:
                simulate(self.system, h, count)
            except Exception as error:
                raise SimulationError(
                    f"Newton failed near t = {done * h:.4f} s with h = {h}. "
                    "Try a smaller time step."
                ) from error
            rec._capture((done + count) * h)
        rec.wall_time = time.perf_counter() - started
        rec.steps = total
        return rec.freeze()


class Recording:
    """A sampled trajectory, ready to animate or plot."""

    def __init__(self, world):
        self.world = world
        self.times = []
        self.positions = []
        self.rotations = []
        self.spins = []
        self.momenta = []
        self.energies = []
        self.kinetic = []
        self.link_ends = []
        self.wall_time = 0.0
        self.steps = 0

    def _capture(self, t):
        bodies = self.world.bodies
        self.times.append(t)
        self.positions.append([b.position for b in bodies])
        self.rotations.append([b.rotation for b in bodies])
        self.spins.append([b.angular_velocity for b in bodies])
        self.momenta.append([b.angular_momentum() for b in bodies])
        self.energies.append(self.world.energy())
        self.kinetic.append(sum(b.kinetic_energy() for b in bodies))
        self.link_ends.append(self.world._link_ends())

    def freeze(self):
        self.times = np.array(self.times)
        self.positions = np.array(self.positions)
        self.rotations = np.array(self.rotations)
        self.spins = np.array(self.spins)
        self.momenta = np.array(self.momenta)
        self.energies = np.array(self.energies)
        self.kinetic = np.array(self.kinetic)
        return self

    @property
    def energy_error(self):
        """Largest energy change, relative to the largest kinetic energy reached."""
        change = np.max(np.abs(self.energies - self.energies[0]))
        scale = max(np.max(self.kinetic), abs(self.energies[0]))
        return change / scale if scale > 0 else change

    def summary(self):
        return (f"{self.steps} steps in {self.wall_time:.2f} s, "
                f"relative energy error {self.energy_error:.2e}")

    def show(self, speed=1.0, **options):
        """Interactive 3D animation; `speed` below 1 plays in slow motion."""
        from .view import show

        return show(self, speed=speed, **options)

    def plot_energy(self, ax=None):
        import matplotlib.pyplot as plt

        ax = ax or plt.subplots(figsize=(7, 2.6))[1]
        ax.plot(self.times, self.energies - self.energies[0], color=PALETTE[1])
        ax.set_xlabel("time [s]")
        ax.set_ylabel(r"$E(t) - E(0)$ [J]")
        ax.set_title("Energy error stays bounded")
        ax.grid(alpha=0.3)
        return ax
