import sys, time, itertools, pickle
sys.path.insert(0, ".")
from oct9d_common import *
t0=time.time()
TW = TW3()
pools = pickle.load(open("pools.pkl","rb"))
poolA, poolB = pools['poolA'], pools['poolB']

def colbasis(M):
    """column basis (sympy) -> list of pivot column vectors"""
    rref, pivs = M.rref()
    return [M[:,i] for i in pivs]

def quotient(Pmat, dim):
    """given product matrix (dim x m), return function class(v) -> quotient coords + rank of P"""
    Bp = colbasis(Pmat)
    r = len(Bp)
    M = sp.Matrix.hstack(*Bp) if Bp else sp.zeros(dim,0)
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
        return x[r:]
    return cls, r, Bp, E

# P for {2,3,5}
pcols = []
for g in Gargs:
    for p in (2,3,5): pcols.append(col_prod(g, p, SS235))
for tr in itertools.combinations_with_replacement((2,3,5),3): pcols.append(col_logprod(tr, SS235))
Pmat235 = sp.Matrix.hstack(*pcols)
cls235, rk235, Bp235, E235 = quotient(Pmat235, 27)
print("P235 dim", rk235, "-> quotient dim", 27-rk235, f"({time.time()-t0:.1f}s)")

tau = {e: cls235(TW[e]) for e in (1,2,3)}
gam = {g: cls235(col_Li3(g, SS235)) for g in poolA}
print("classes computed", f"({time.time()-t0:.1f}s)")

# Greedy: cover T = span(tau1,tau2,tau3) with few candidate classes
T = sp.Matrix.hstack(tau[1], tau[2], tau[3])
rkT = T.rank()
print("rank T =", rkT)
def dimcap(S):
    if not S: return 0
    MS = sp.Matrix.hstack(*S)
    return rkT + MS.rank() - sp.Matrix.hstack(MS, T).rank()
chosen = []
cand = list(poolA)
# prefer edge points and small heights: sort candidates by (height)
def height(z): return max(abs(z).numerator, abs(z).denominator)
cand.sort(key=lambda z: (height(z), str(z)))
while dimcap(chosen) < rkT:
    best, bestd = None, dimcap(chosen)
    for g in cand:
        if g in [c[0] for c in chosen]: continue
        d = dimcap([c[1] for c in chosen] + [gam[g]])
        if d > bestd:
            # tie-break by height handled by iteration order
            if best is None or d > bestd:
                best, bestd = (g, gam[g]), d
        if d == rkT: 
            best, bestd = (g, gam[g]), d
            break
    if best is None:
        print("STUCK at dimcap", dimcap(chosen)); break
    chosen.append(best)
    print(f"add {best[0]} (height {height(best[0])}): dimcap -> {bestd}  ({time.time()-t0:.1f}s)")
print("chosen for edges 1-3:", [c[0] for c in chosen])
pickle.dump({'chosen123': [c[0] for c in chosen]}, open("short123.pkl","wb"))
