# ASC-ODE-EIGEN
Rigid Body Simulation Engine implementation with own BLA module.
A python notebook showing how to use the library will come soon.
This simulation engine is part of my Bachelor Thesis which can also be found in this repository.
A JupyterLite notebook can be found at https://linux-03.github.io/rigid_body_inter/. Please note that this website is still under development, and the system is not yet running at full efficiency. As a result, the code may run more slowly and additional example content may be displayed.


## Conventions

- `p_trans` is the world-frame linear momentum, `m * v`.
- `p_skew` is the body-frame momentum, `1/2 * J * omega` (thesis eq. 4.13).
- Gravity defaults to `(0, -9.81, 0)` and is scaled by each body's mass. Set it with `RigidBodySystem::Gravity()` in C++ or `rbs.gravity` in Python.
- Add a body to the system before creating connectors that refer to it.

## Tests

```
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```
