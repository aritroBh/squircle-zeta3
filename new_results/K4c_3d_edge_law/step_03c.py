"""step_03c.py -- exact dihedral cosines for model (b), outward normals."""
import sympy as sp
from sympy import sqrt, simplify, acos, pi, nsimplify

a = (4*sqrt(3)-3)/13
c = (2*sqrt(3)-3)/2
A = (a,0,0); P = (c,c,0); P3 = (c,0,c); P4 = (c,0,-c); R = (0,c,c)
def sub(u,v): return tuple(u[i]-v[i] for i in range(3))
def cross(u,v): return (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
def dot(u,v): return sp.simplify(sum(u[i]*v[i] for i in range(3)))
def norm2(u): return sp.simplify(dot(u,u))

# orbit1: edge A-P, faces (A,P,P3), (A,P,P4); outward normals
nL = cross(sub(P,A), sub(P3,A))   # (c^2, c(a-c), c(a-c)) outward
nR = cross(sub(P4,A), sub(P,A))
cosang1 = sp.simplify(dot(nL,nR)/sp.sqrt(norm2(nL)*norm2(nR)))
print("orbit1: cos(angle between outward normals) =", sp.simplify(cosang1))
print("        = numeric", sp.N(cosang1,30))
print("        cos(dihedral) =", sp.simplify(-cosang1), "numeric", sp.N(-cosang1,30))
print("        dihedral =", sp.N(pi-acos(cosang1),30))
print("        minimal poly of cos(dih):", sp.minpoly(sp.N(-cosang1,60)) if False else "")

# orbit2: edge P-P3, faces (A,P,P3), (R,P,P3)
nT = cross(sub(P3,P), sub(R,P))
nT = tuple(-x for x in nT)  # outward (+,+,+)
cosang2 = sp.simplify(dot(nL,nT)/sp.sqrt(norm2(nL)*norm2(nT)))
print("orbit2: cos(angle between outward normals) =", sp.simplify(cosang2))
print("        cos(dihedral) =", sp.simplify(-cosang2), "numeric", sp.N(-cosang2,30))
print("        dihedral =", sp.N(pi-acos(cosang2),30))
# tan of half-dihedral maybe cleaner
import mpmath as mp
mp.mp.dps = 50
d1 = mp.pi - mp.acos(mpf(sp.N(cosang1,50)))
d2 = mp.pi - mp.acos(mpf(sp.N(cosang2,50)))
print("orbit1 dihedral 50d:", mp.nstr(d1,50))
print("orbit2 dihedral 50d:", mp.nstr(d2,50))
print("sin/cos half angles:")
print("  orbit1 cos(dih) 50d:", mp.nstr(mp.cos(d1),50), " sin:", mp.nstr(mp.sin(d1),50))
print("  orbit2 cos(dih) 50d:", mp.nstr(mp.cos(d2),50), " sin:", mp.nstr(mp.sin(d2),50))
# nsimplify the cosines on Q(sqrt3)
c1 = nsimplify(mpf(float(mp.cos(d1))), [sqrt(3)])
c2 = nsimplify(mpf(float(mp.cos(d2))), [sqrt(3)])
print("cos(dih1) in Q(sqrt3)?", c1)
print("cos(dih2) in Q(sqrt3)?", c2)
