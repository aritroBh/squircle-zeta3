import sys, time, itertools, pickle
sys.path.insert(0, ".")
from oct9d_common import *
t0=time.time()
TW = TW3()
pools = pickle.load(open("pools.pkl","rb"))
poolA, poolB = pools['poolA'], pools['poolB']

def colbasis(M):
    rref, pivs = M.rref()
    return [M[:,i] for i in pivs]

def quotient(Pmat, dim):
    Bp = colbasis(Pmat)
    r = len(Bp)
    M = sp.Matrix.hstack(*Bp)
    E = []
    for i in range(dim):
        ei = sp.zeros(dim,1); ei[i,0]=1
        if (sp.Matrix.hstack(M, *E, ei)).rank() > r + len(E):
            E.append(ei)
        if r + len(E) == dim: break
    Mfull = sp.Matrix.hstack(M, *E)
    Minv = Mfull.inv()
    def cls(v):
        x = Minv*sp.Matrix(v)
        return x[r:, 0]
    return cls, r

pcols = []
for g in Gargs:
    for p in (2,3,5): pcols.append(col_prod(g, p, SS235))
for tr in itertools.combinations_with_replacement((2,3,5),3): pcols.append(col_logprod(tr, SS235))
Pmat235 = sp.Matrix.hstack(*pcols)
cls235, rk235 = quotient(Pmat235, 27)
tau = {e: cls235(TW[e]) for e in (1,2,3)}
gam = {g: cls235(col_Li3(g, SS235)) for g in poolA}
T = sp.Matrix.hstack(tau[1], tau[2], tau[3]); rkT = T.rank()
def dimcap(S):
    if not S: return 0
    MS = sp.Matrix.hstack(*S)
    return rkT + MS.rank() - sp.Matrix.hstack(MS, T).rank()
def height(z): return max(abs(z).numerator, abs(z).denominator)
cand = sorted(poolA, key=lambda z: (height(z), str(z)))
chosen = []
while dimcap([c[1] for c in chosen]) < rkT:
    cur = dimcap([c[1] for c in chosen])
    best = None
    for g in cand:
        if any(g == c[0] for c in chosen): continue
        d = dimcap([c[1] for c in chosen] + [gam[g]])
        if d > cur:
            best = (g, gam[g], d)
            if d == rkT: break
    if best is None: print("STUCK at", cur); break
    chosen.append((best[0], best[1]))
    print(f"add z={best[0]} (h={height(best[0])}): dimcap -> {best[2]}")
print("chosen for edges 1-3 (minimal %d):"%len(chosen), [c[0] for c in chosen])
pickle.dump({'chosen123': [c[0] for c in chosen]}, open("short123.pkl","wb"))
print(f"({time.time()-t0:.1f}s)")
