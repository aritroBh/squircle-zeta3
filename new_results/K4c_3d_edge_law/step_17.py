
step03b = r'''"""step_03b.py -- exact sympy geometry for model (b): faces, areas, dihedrals, S(1)."""
import sympy as sp
from sympy import sqrt, Rational, simplify, nsimplify, acos, pi

a = (4*sqrt(3)-3)/13       # axis vertex distance
c = (2*sqrt(3)-3)/2        # diag vertex coord
# vertices
V = []
for i in range(3):
    for s in [1,-1]:
        v = [0,0,0]; v[i] = s*a; V.append(tuple(v))
for perm in [(0,1,2),(1,2,0),(2,0,1)]:
    for s1 in [1,-1]:
        for s2 in [1,-1]:
            v = [0,0,0]; v[perm[1]]=s1*c; v[perm[2]]=s2*c; V.append(tuple(v))
def sub(u,v): return tuple(u[i]-v[i] for i in range(3))
def cross(u,v): return (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
def dot(u,v): return sum(u[i]*v[i] for i in range(3))
# edge orbits (from step_03): orbit1 axis-diag L1^2=(a-c)^2+c^2; orbit2 diag-diag L2 = c*sqrt2
L1sq = sp.simplify((a-c)**2 + c**2)
L2sq = sp.simplify(2*c**2)
print("orbit1 (axis-diag) L^2 =", sp.simplify(L1sq), "=", sp.nsimplify(sp.N(L1sq,30)), " numeric:", sp.N(L1sq,20))
print("orbit2 (diag-diag) L^2 =", sp.simplify(L2sq), "  L2 =", sp.simplify(sp.sqrt(L2sq)), " numeric:", sp.N(sp.sqrt(L2sq),20))
print("orbit1 L1 =", sp.nsimplify(sp.N(sp.sqrt(L1sq),40)), " numeric:", sp.N(sp.sqrt(L1sq),20))
# face type 1: (a,0,0),(c,c,0),(c,0,c)
A = (a,0,0); P = (c,c,0); Q = (c,0,c)
n1 = cross(sub(P,A), sub(Q,A))
A1sq = sp.simplify(dot(n1,n1)/4)
print("type-1 face area^2 =", sp.simplify(A1sq), " numeric:", sp.N(A1sq,20))
A1 = sp.sqrt(A1sq)
print("type-1 face area =", sp.simplify(A1), " numeric:", sp.N(A1,20))
# face type 2: equilateral side L2
A2 = sp.sqrt(3)/4*L2sq
print("type-2 face area =", sp.simplify(A2), " numeric:", sp.N(A2,20))
S = sp.simplify(24*A1 + 8*A2)
print("S(1) =", sp.simplify(S))
print("S(1) numeric:", sp.N(S,40))
# dihedral angle orbit1: edge (a,0,0)-(c,c,0): faces (A,P,Q) and (A,P,(c,-c,0))
Q2 = (c,-c,0)
n1b = cross(sub(P,A), sub(Q2,A))
cosd1 = sp.simplify(-dot(n1,n1b)/sp.sqrt(dot(n1,n1)*dot(n1b,n1b)))
print("orbit1 cos(dihedral) =", sp.simplify(cosd1), " numeric:", sp.N(cosd1,20))
print("       dihedral =", sp.N(acos(-cosd1* -1 if False else cosd1),5))
dih1 = acos(-cosd1) if cosd1.evalf()<0 else None
# interior dihedral = pi - angle between outward normals; outward normals n1,n1b already outward? check sign
# just print both
print("       angle between normals:", sp.N(acos(sp.simplify(dot(n1,n1b)/sp.sqrt(dot(n1,n1)*dot(n1b,n1b)))),20))
print("       dihedral (pi - that):", sp.N(pi - acos(sp.simplify(dot(n1,n1b)/sp.sqrt(dot(n1,n1)*dot(n1b,n1b)))),20))
# orbit2: edge (c,c,0)-(c,0,c): faces (A,P,Q) and ((0,c,c),P,Q)
R = (0,c,c)
n2b = cross(sub(Q,P), sub(R,P))
n2a = cross(sub(Q,P), sub(A,P))
c22 = sp.simplify(dot(n2a,n2a)); c21 = sp.simplify(dot(n2a,n2b)); c11 = sp.simplify(dot(n2b,n2b))
ang = sp.simplify(c21/sp.sqrt(c22*c11))
print("orbit2 cos(angle between normals) =", sp.simplify(ang), " numeric:", sp.N(ang,20))
print("       dihedral (pi - that):", sp.N(pi - acos(ang),20))
# check outward orientation sign: face (A,P,Q): centroid ~ (a/3+2c/3, ...)>0, n1 should point outward (positive dot)
cen = tuple((A[i]+P[i]+Q[i])/3 for i in range(3))
print("outward check n1.centroid =", sp.N(dot(n1,cen),5), " n2a.centroid(face A,P,Q) =", sp.N(dot(n2a,cen),5))
'''
with open('step_03b.py','w') as f:
    f.write(step03b)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_03b.py'],
                   capture_output=True, text=True)
print(r.stdout)
print(r.stderr[-3000:] if r.returncode else "exit 0")
