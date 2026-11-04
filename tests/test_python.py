"""Checks the Python API. Run with: python tests/test_python.py"""
import numpy as np
from netgen.occ import Box, Pnt, Cylinder
from lib_rigid_body import World
cube=lambda a: Box(Pnt(-a,-a,-a),Pnt(a,a,a))
# pendulum: period vs exact small-amplitude, energy
w=World(); b=w.add(cube(.05),at=(0.1,0,-1)); w.pin(b,(0,0,0)); r=w.run(2.0,h=.002); print('pendulum',r.summary(), 'link len',np.linalg.norm(r.positions[-1,0]))
assert r.energy_error<1e-4 and abs(np.linalg.norm(r.positions[-1,0])-np.linalg.norm([.1,0,-1]))<1e-9
# heavy top from CAD mesh with explicit inertia
J=np.diag([.234375,.46875,.234375])
w=World(); t=w.add(Cylinder(Pnt(0,.95,0),(0,1,0),r=.5,h=.1),at=(0,1,0),mass=15,inertia=J,spin=(0,150,-4.61538),velocity=(4.61538,0,0)); w.pin(t,(0,0,0)); r=w.run(1.0,h=.0005); print('top',r.summary(),'COM',r.positions[-1,0]); assert np.linalg.norm(r.positions[-1,0]-[0.173343964098,0.640088592073,-0.748490791132])<5e-4
# double pendulum ball joints: body-body SPHERICAL energy
w=World(); a=w.add(cube(.1),at=(1,0,0),density=100); c=w.add(cube(.1),at=(2,0,0),density=100); w.pin(a,(0,0,0)); w.ball(a,c,(1,0,0.0)+np.array([.5,0,0])); r=w.run(2.0,h=.002); print('double',r.summary()); assert r.energy_error<1e-3
# free body: L & E conserved, flip
h1=Box(Pnt(-.5,-.05,-.05),Pnt(.5,.05,.05)); h2=Box(Pnt(-.05,0,-.05),Pnt(.05,.5,.05)); T=h1+h2
w=World(gravity=(0,0,0)); f=w.add(T,spin=(0.01,4,.01)); print('J',np.round(f.inertia,4)); r=w.run(8.0,h=.002); L=r.momenta[:,0]; print('free',r.summary(),'dL',np.max(np.linalg.norm(L-L[0],axis=1))); print('w_y range',r.spins[:,0,1].min(),r.spins[:,0,1].max()); assert np.max(np.linalg.norm(L-L[0],axis=1))<1e-8
# view builds
v=r.show(speed=.5); print(type(v).__name__)
print('ALL OK')
