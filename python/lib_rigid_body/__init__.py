"""Structure-preserving rigid body simulation with a small, friendly Python API.

    from lib_rigid_body import World
    from netgen.occ import Box, Pnt

    world = World()
    cube = world.add(Box(Pnt(-.2, -.2, -.2), Pnt(.2, .2, .2)), at=(1, 0, 0))
    world.pin(cube, (0, 0, 0))
    world.run(duration=4).show()
"""
from . import rigid_body_FEM as core
from .rigid_body_FEM import bla
from .world import Body, Recording, SimulationError, World

__all__ = ["Body", "Recording", "SimulationError", "World", "bla", "core"]
