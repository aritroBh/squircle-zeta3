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
Pmat235 = sp.Matrix.hstack(*pcols)
cls235, rkP235 = quotient(Pmat235, 27)
tau = {e: cls235(TW[e]) for e in (1,2,3)}

# niceness classification
EP = set()
for e in (1,2,3):
    for z in edge_points(e): EP.add(z)
ORB = set()
for z in list(EP): ORB |= orbit(z)
SM = set(SMALL235)
ORBSM = set()
for z in list(SM): ORBSM |= orbit(z)
def nice(g):
    g = Fr(g)
    if g in EP: return 0
    if g in ORB: return 1
    if g in SM: return 2
    if g in ORBSM: return 3
    return 4

gamQ = {g: cls235(col_Li3(sp.Rational(g), SS235)) for g in Gargs}
gamQ = {g: v for g, v in gamQ.items() if any(v)}
keys = sorted(gamQ.keys(), key=lambda g: (nice(g), max(abs(Fr(g).numerator), abs(Fr(g).denominator)), str(g)))
print("candidates by niceness:", [(str(g), nice(g)) for g in keys[:25]])

M8 = {g: sp.Matrix(gamQ[g]) for g in keys}
def det8(sel):
    return abs(sp.Matrix.hstack(*[M8[g] for g in sel]).det())

# greedy build from nicest, then swap-optimize det
random.seed(1)
best = None
def try_improve(sel, d):
    improved = True
    while improved:
        improved = False
        for i in range(len(sel)):
            for g in keys:
                if g in sel: continue
                sel2 = sel[:i] + [g] + sel[i+1:]
                try: d2 = det8(sel2)
                except Exception: continue
                if d2 != 0 and d2 < d:
                    sel, d = sel2, d2; improved = True
    return sel, d
# starts: greedy by niceness with randomization
for trial in range(60):
    if trial == 0:
        sel = []
        for g in keys:
            if len(sel) == 8: break
            if sp.Matrix.hstack(*[M8[x] for x in sel]+[M8[g]]).rank() > len(sel): sel.append(g)
    else:
        pool0 = [g for g in keys if nice(g) <= 2]
        sel = []
        cand = pool0 + [g for g in keys if nice(g) == 3]
        while len(sel) < 8:
            g = random.choice(cand[:min(len(cand), 30)])
            if g in sel: continue
            if sp.Matrix.hstack(*[M8[x] for x in sel]+[M8[g]]).rank() > len(sel): sel.append(g)
    d = det8(sel)
    sel, d = try_improve(sel, d)
    if best is None or d < best[1]:
        best = (sel, d)
        print(f"trial {trial}: det {d} = {sp.nsimplify(d)} <- { [str(g) for g in sel] }  ({time.time()-t0:.1f}s)")
print("BEST:", best[1], [str(g) for g in best[0]])
pickle.dump({'basis235': best[0], 'det': best[1]}, open("basis235.pkl","wb"))
