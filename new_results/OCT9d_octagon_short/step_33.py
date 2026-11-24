import sys, time, itertools, pickle, random
import numpy as np
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
for g in Gargs:
    for p in (2,3,5): pcols.append(col_prod(g, p, SS235))
for tr in itertools.combinations_with_replacement((2,3,5),3): pcols.append(col_logprod(tr, SS235))
cls235, rkP235 = quotient(sp.Matrix.hstack(*pcols), 27)
tau = {e: cls235(TW[e]) for e in (1,2,3)}
EP = set()
for e in (1,2,3):
    for z in edge_points(e): EP.add(z)
ORB = set()
for z in list(EP): ORB |= orbit(z)
ORBSM = set()
for z in list(SMALL235): ORBSM |= orbit(z)
def nice(g):
    g = Fr(g)
    if g in EP: return 0
    if g in ORB: return 1
    if g in SMALL235: return 2
    if g in ORBSM: return 3
    return 4
gamQ = {g: cls235(col_Li3(sp.Rational(g), SS235)) for g in Gargs}
gamQ = {g: v for g, v in gamQ.items() if any(v)}
keys = sorted(gamQ.keys(), key=lambda g: (nice(g), max(abs(Fr(g).numerator), abs(Fr(g).denominator)), str(g)))
M8 = {g: sp.Matrix(gamQ[g]) for g in keys}
def det8(sel): return abs(sp.Matrix.hstack(*[M8[g] for g in sel]).det())
random.seed(2)
best = (None, None)
for trial in range(400):
    if trial == 0:
        sel = []
        for g in keys:
            if len(sel) == 8: break
            if sp.Matrix.hstack(*[M8[x] for x in sel]+[M8[g]]).rank() > len(sel): sel.append(g)
    else:
        sel = random.sample(keys, 8)
        if sp.Matrix.hstack(*[M8[x] for x in sel]).rank() < 8: continue
    d = det8(sel)
    improved = True
    while improved:
        improved = False
        for i in range(8):
            for g in keys:
                if g in sel: continue
                sel2 = sel[:i]+[g]+sel[i+1:]
                d2 = det8(sel2)
                if d2 != 0 and d2 < d: sel, d = sel2, d2; improved = True
    if best[1] is None or d < best[1]:
        best = (list(sel), d)
        print(f"trial {trial}: det {d} <- {[str(g) for g in sel]} ({time.time()-t0:.1f}s)")
    if best[1] == 1: break
print("BEST det:", best[1], [str(g) for g in best[0]])
pickle.dump({'basis235': best[0], 'det': best[1]}, open("basis235.pkl","wb"))
