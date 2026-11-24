import sys, time, itertools, pickle, random
sys.path.insert(0, ".")
from oct9d_common import *
t0=time.time()
TW = TW3()
def colbasis(M):
    rref, pivs = M.rref(); return [M[:,i] for i in pivs]
def quotient(Pmat, dim):
    Bp = colbasis(Pmat); r = len(Bp)
    M = sp.Matrix.hstack(*Bp); E = []
    for i in range(dim):
        ei = sp.zeros(dim,1); ei[i,0]=1
        if (sp.Matrix.hstack(M, *E, ei)).rank() > r+len(E): E.append(ei)
        if r+len(E) == dim: break
    Minv = (sp.Matrix.hstack(M, *E)).inv()
    def cls(v):
        x = Minv*sp.Matrix(v); return x[r:,0]
    return cls, r
pcols = []
for w in negw7:
    for p in (2,3,5,7): pcols.append(col_prod(sp.Rational(w), p, SS2357))
for tr in itertools.combinations_with_replacement((2,3,5,7),3): pcols.append(col_logprod(tr, SS2357))
Pmat0 = sp.Matrix.hstack(*pcols)
cls0, rkP0 = quotient(Pmat0, 64)

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

gamB = {}
for g in sols7:
    v = cls0(col_Li3(sp.Rational(g), SS2357))
    if any(v): gamB[g] = v
keys = sorted(gamB.keys(), key=lambda g:(nice0(g), max(abs(Fr(g).numerator),abs(Fr(g).denominator)), str(g)))
print("nicest 30:", [(str(g), nice0(g)) for g in keys[:30]])
M20 = {g: sp.Matrix(gamB[g]) for g in keys}
def det20(sel): return abs(sp.Matrix.hstack(*[M20[g] for g in sel]).det())
sel = []
for g in keys:
    if len(sel) == 20: break
    if sp.Matrix.hstack(*[M20[x] for x in sel]+[M20[g]]).rank() > len(sel): sel.append(g)
d = det20(sel)
print("greedy det:", d, f"({time.time()-t0:.1f}s)")
swc = keys[:130]
improved = True; it = 0
while improved and it < 3:
    improved = False; it += 1
    for i in range(20):
        for g in swc:
            if g in sel: continue
            sel2 = sel[:i]+[g]+sel[i+1:]
            d2 = det20(sel2)
            if d2 != 0 and d2 < d: sel, d = sel2, d2; improved = True
        print(f"pass {it} i={i}: det {d} ({time.time()-t0:.1f}s)")
        if d == 2592: break
    if d == 2592: break
print("FINAL det:", d)
print("basis:", [str(g) for g in sel])
print("niceness:", [nice0(g) for g in sel])
pickle.dump({'basis0': [str(g) for g in sel], 'det': str(d)}, open("basis0.pkl","wb"))
