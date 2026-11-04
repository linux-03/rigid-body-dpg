# Rigid body DPG

A structure-preserving rigid body simulator, written for my bachelor thesis ([`bachelor.pdf`](bachelor.pdf)). It implements a second-order symplectic time integrator on the Lie group R^3 x SO(3), derived from a discontinuous Petrov-Galerkin discretization of the Hamilton-Pontryagin principle. The linear algebra and automatic differentiation are self-written.

## Try it online

**[Open the examples](https://linux-03.github.io/rigid-body-dpg/lab/index.html?path=01_pendulum.ipynb)**. They run in your browser, with nothing to install.

```python
from netgen.occ import Box, Pnt
from lib_rigid_body import World

world = World()
cube = world.add(Box(Pnt(-.25, -.25, -.25), Pnt(.25, .25, .25)), at=(1.5, 0, 0), density=100)
world.pin(cube, (0, 0, 0))
world.run(duration=6, h=0.005).show()
```

To run the examples locally:

```
git -c url."https://github.com/".insteadOf="git@github.com:" submodule update --init --recursive
pip install . netgen-mesher pythreejs jupyterlab
jupyter lab examples
```

## Conventions

- `p_trans` is the world-frame linear momentum, `m * v`.
- `p_skew` is the body-frame momentum, `1/2 * J * omega` (thesis eq. 4.13).
- Gravity is scaled by each body's mass. Python's `World` defaults to `(0, 0, -9.81)`. The C++ `RigidBodySystem` defaults to `(0, -9.81, 0)`; set it with `Gravity()`.
- Add a body to the system before creating connectors that refer to it.

## Tests

```
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```
