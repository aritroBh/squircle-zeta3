
step03 = r'''"""step_03.py -- exact geometry of K_1 for models (a) and (b)."""
import os, itertools
os.chdir("research/08_theory/K4c_3d_edge_law")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 60
import numpy as np
import sympy as sp
from sympy import sqrt, Rational, acos, simplify, nsimplify

sq2, sq3 = mp.sqrt(2), mp.sqrt(3)

# ---------------- MODEL (a): tetrakis hexahedron ----------------
print("="*70); print("MODEL (a) exact"); print("="*70)
pA = 1/(2*sq2)   # axis vertices
pB = 1/(3*sq2)   # body vertices
VA = []
for i in range(3):
    for s in [1,-1]:
        v = [mpf(0)]*3; v[i] = s*pA; VA.append(tuple(v))
for s1,s2,s3 in itertools.product([1,-1],repeat=3):
    VA.append((s1*pB, s2*pB, s3*pB))
VA = [np.array(v,float) for v in VA]
num = np.load("K1V_A.npy")
# match
ok = all(min(np.linalg.norm(v-w) for w in num) < 1e-8 for v in VA)
print("all 14 exact vertices match numeric hull:", ok)
# faces: triangles; from earlier: facet from Z vertex vZ=(sq2,0,2sq2): vertices (0,0,pA),(pB,pB,pB),(pB,-pB,pB)
# face = {axis apex + two adjacent body verts} and {??}: count via numeric adjacency from step_02
# edges: orbit1: apex-body L=1/(2sqrt2)=0.353553; orbit2: body-body L=sqrt2/3=0.471405
L1 = float(1/(2*sq2)); L2 = float(sq2/3)
def d(v,w): return np.linalg.norm(np.array(v)-np.array(w))
n1 = sum(1 for i in range(14) for j in range(i+1,14) if abs(d(VA[i],VA[j])-L1)<1e-9)
n2 = sum(1 for i in range(14) for j in range(i+1,14) if abs(d(VA[i],VA[j])-L2)<1e-9)
print(f"edges with L=1/(2sqrt2): {n1}, with L=sqrt2/3: {n2}")
# dihedral: face normals; faces: apex(i,s) connects to 4 body verts forming a square cycle
# take faces: (e3 apex, (pB,pB,pB),(pB,-pB,pB)) and (e3 apex,(pB,-pB,pB),(-pB,-pB,pB)): share edge apex-(pB,-pB,pB) [orbit1]
A1=np.array([0,0,pA]); P1=np.array([pB,pB,pB]); Q1=np.array([pB,-pB,pB]); R1=np.array([-pB,-pB,pB])
def fnrm(a,b,c):
    n=np.cross(b-a,c-a); return n/np.linalg.norm(n)
nA=fnrm(A1,P1,Q1); nB=fnrm(A1,Q1,R1)
dih1 = np.pi - np.arccos(np.dot(nA,nB))
print(f"orbit1 (apex-body) dihedral: {dih1:.15f}  arccos(-4/5)={float(mp.acos(mpf(-4)/5)):.15f}")
# orbit2 edge: body-body e.g. (pB,pB,pB)-(pB,-pB,pB): shared by faces (e3apex,P1,Q1) and (?? other apex?)
# body edge P1-Q1: adjacent faces have apexes e3 (z+) and e1? (x+): face (apex_x, P1, Q1) with apex_x=(pA,0,0)
A2=np.array([pA,0,0])
nC=fnrm(A2,P1,Q1)
dih2 = np.pi - np.arccos(np.dot(nA,nC))
print(f"orbit2 (body-body) dihedral: {dih2:.15f}  arccos(-4/5)={float(mp.acos(mpf(-4)/5)):.15f}")
# exact area check
import sympy.geometry as sg
# triangle (A1,P1,Q1) exact
a1=sp.Point3D(0,0,1/(2*sqrt(2))); p1=sp.Point3D(1/(3*sqrt(2)),1/(3*sqrt(2)),1/(3*sqrt(2)))
q1=sp.Point3D(1/(3*sqrt(2)),-1/(3*sqrt(2)),1/(3*sqrt(2)))
ar = sg.Triangle(a1,p1,q1).area
print("one face area exact:", sp.simplify(ar), "= S/24? S =", sp.simplify(24*ar), " 2sqrt5/3 =", 2*sqrt(5)/3, "equal:", sp.simplify(24*ar-2*sqrt(5)/3)==0)

# ---------------- MODEL (b) ----------------
print("="*70); print("MODEL (b) exact"); print("="*70)
av = (4*sq3-3)/13    # axis vertex distance
cv = (2*sq3-3)/2     # diag vertex coord
VB = []
for i in range(3):
    for s in [1,-1]:
        v=[mpf(0)]*3; v[i]=s*av; VB.append(tuple(v))
for perm in [(0,1,2),(1,2,0),(2,0,1)]:  # (0,+-c,+-c),(+-c,0,+-c),(+-c,+-c,0)
    for s1,s2 in itertools.product([1,-1],repeat=2):
        v=[mpf(0)]*3
        v[perm[1]]=s1*cv; v[perm[2]]=s2*cv
        VB.append(tuple(v))
VB=[np.array(v,float) for v in VB]
numB = np.load("K1V_B.npy")
okB = all(min(np.linalg.norm(v-w) for w in numB) < 1e-8 for v in VB) and \
      all(min(np.linalg.norm(v-w) for v in VB) < 1e-8 for w in numB)
print("all 18 exact vertices match numeric hull:", okB)
# edges: find pairs at the two orbit lengths
LA, LB_ = 0.242413295, 0.328169399
pairs = {}
for i in range(18):
    for j in range(i+1,18):
        dd = d(VB[i],VB[j])
        if abs(dd-LA)<1e-6: pairs[(i,j)]=('orbit1',dd)
        elif abs(dd-LB_)<1e-6: pairs[(i,j)]=('orbit2',dd)
print(f"found {len(pairs)} edges")
from collections import Counter
def vtype(i): return 'axis' if i<6 else 'diag'
cnt = Counter((v[0], tuple(sorted((vtype(k[0]),vtype(k[1]))))) for k,v in pairs.items())
print("edge orbits by vertex types:", dict(cnt))
# exact lengths
mp.mp.dps = 60
ex1 = nsimplify(LA, [sqrt(2), sqrt(3), sqrt(6)])
ex2 = nsimplify(LB_, [sqrt(2), sqrt(3), sqrt(6)])
print("orbit1 L =", ex1, "=", sp.N(ex1,20))
print("orbit2 L =", ex2, "=", sp.N(ex2,20))
# refine lengths at 60 dps then nsimplify
def d60(v,w): return mp.sqrt(sum((mpf(v[i])-mpf(w[i]))**2 for i in range(3)))
s1 = sorted(set(k for k,v in pairs.items() if v[0]=='orbit1'))
s2 = sorted(set(k for k,v in pairs.items() if v[0]=='orbit2'))
L1h = d60(VB[s1[0][0]], VB[s1[0][1]]); L2h = d60(VB[s2[0][0]], VB[s2[0][1]])
print("orbit1 L 60d:", mp.nstr(L1h,50), "->", nsimplify(mpf(float(L1h)),[sqrt(2),sqrt(3),sqrt(6)]))
print("orbit2 L 60d:", mp.nstr(L2h,50), "->", nsimplify(mpf(float(L2h)),[sqrt(2),sqrt(3),sqrt(6)]))
# dihedrals: need faces: a face = 3 vertices lying on plane <y,vZ>=1. Reconstruct faces from adjacency:
# face = triple of pairwise-adjacent vertices (all pairs among orbit edges)
adj = {}
for (i,j) in pairs:
    adj.setdefault(i,set()).add(j); adj.setdefault(j,set()).add(i)
faces = set()
for i in range(18):
    for j,k in itertools.combinations(sorted(adj.get(i,[])),2):
        if k in adj.get(j,set()):
            faces.add(tuple(sorted((i,j,k))))
faces = sorted(faces)
print(f"{len(faces)} triangular faces")
vf = {}
for fi,f in enumerate(faces):
    for vi in f: vf.setdefault(vi,[]).append(fi)
val = Counter(len(vf[i]) for i in range(18))
print("vertex valences:", dict(val))
fn = [fnrm(VB[f[0]],VB[f[1]],VB[f[2]]) for f in faces]
dihs = []
for (i,j),(orb,_) in pairs.items():
    sh = [fi for fi in vf[i] if fi in vf[j]]
    n1_,n2_ = fn[sh[0]], fn[sh[1]]
    dd = np.pi - np.arccos(np.clip(np.dot(n1_,n2_),-1,1))
    dihs.append((orb, dd))
for orb in ['orbit1','orbit2']:
    vals = [x for o,x in dihs if o==orb]
    print(orb, "dihedral:", f"{vals[0]:.15f}", "count", len(vals), "cos=", f"{np.cos(vals[0]):.15f}",
          "nsimplify cos:", nsimplify(np.cos(vals[0]), [sqrt(2),sqrt(3),sqrt(6),sqrt(5)]))
# surface area exact
SB = mpf(0)
for f in faces:
    a,b,c = [VB[i] for i in f]
    cr = np.cross(b-a,c-a)
    SB += mpf(0.5)*mpf(float(np.linalg.norm(cr)))
print("S(1) model b 60d:", mp.nstr(SB,50))
print("nsimplify:", nsimplify(mpf(float(SB)), [sqrt(2), sqrt(3), sqrt(6)], full=True))
'''
with open('step_03.py','w') as f:
    f.write(step03)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_03.py'],
                   capture_output=True, text=True)
print(r.stdout)
print(r.stderr[-3000:] if r.returncode else "exit 0")
