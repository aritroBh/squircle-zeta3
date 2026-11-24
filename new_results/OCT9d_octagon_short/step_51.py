import sys, time, itertools, pickle
import mpmath as mp
sys.path.insert(0, ".")
from oct9d_common import *
from oct9c_common import atom_values
t0=time.time()
TW = TW3()
basis0 = [sp.Rational(g) for g in pickle.load(open("basis0.pkl","rb"))['basis0']]
Gsym0 = {g: col_Li3(sp.Rational(g), SS2357) for g in sols7}
G20 = sp.Matrix.hstack(*[Gsym0[g] for g in basis0])

EP0 = set()
for z in edge_points(0): EP0.add(z)
EP123 = set()
for e in (1,2,3):
    for z in edge_points(e): EP123.add(z)
ORB0 = set()
for z in list(EP0): ORB0 |= orbit(z)
ORB123 = set()
for z in list(EP123): ORB123 |= orbit(z)
ORBSM = set()
for z in list(SMALL235+SMALL7): ORBSM |= orbit(z)
def nice0(g):
    g = Fr(g)
    if g in EP0: return 0
    if g in ORB0: return 1
    if g in EP123 or g in ORB123: return 2
    if g in SMALL235 or g in SMALL7: return 3
    if g in ORBSM: return 4
    return 5
worder = sorted(set([sp.Rational(w) for w in negw7]), key=lambda g:(nice0(g), max(abs(Fr(g).numerator),abs(Fr(g).denominator)), str(g)))
PB0 = []
cur = sp.Matrix.hstack(G20, sp.zeros(64,0))
for w in worder:
    for p in (2,3,5,7):
        c = col_prod(w, p, SS2357)
        if sp.Matrix.hstack(cur, c).rank() > cur.shape[1]:
            PB0.append((w,p)); cur = sp.Matrix.hstack(cur, c)
        if cur.shape[1] == 64: break
    if cur.shape[1] == 64: break
print("PB0 count:", len(PB0), f"({time.time()-t0:.1f}s)")
print("PB0:", [(str(w),p) for w,p in PB0])
M = sp.Matrix.hstack(G20, *[col_prod(w,p,SS2357) for w,p in PB0])
print("det:", M.det(), f"({time.time()-t0:.1f}s)")
x = M.LUsolve(sp.Matrix(TW[0]))
ok = M*x == sp.Matrix(TW[0])
print("solve exact (M x = t_0):", ok, f"({time.time()-t0:.1f}s)")
cl = [sp.Rational(c) for c in x[:20]]; cp = [sp.Rational(c) for c in x[20:]]
print("Li3 coeffs:", [str(c) for c in cl])
print("prod coeffs:", [str(c) for c in cp])
pickle.dump({'PB0': [(str(w),p) for w,p in PB0], 'sol0': [str(c) for c in x]}, open("sol0_short.pkl","wb"))
