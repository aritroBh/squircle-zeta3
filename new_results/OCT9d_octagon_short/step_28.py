import sys, time, itertools, pickle
import numpy as np
sys.path.insert(0, ".")
from oct9d_common import *
t0=time.time()
TW = TW3()
pools = pickle.load(open("pools.pkl","rb"))
poolA, poolB = pools['poolA'], pools['poolB']

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

def intvec(v):
    den = sp.ilcm(*[x.q for x in v])
    return tuple(int(sp.Rational(x)*den) for x in v)
P = 2**31-1
def grank_rows(rows):  # rows: list of int-tuples
    M = np.array([[x % P for x in v] for v in rows], dtype=np.int64) % P
    n, m = M.shape; r = 0
    for c in range(m):
        piv = -1
        for i in range(r, n):
            if M[i,c]: piv = i; break
        if piv < 0: continue
        if piv != r: M[[r,piv]] = M[[piv,r]]
        M[r] = (M[r]*pow(int(M[r,c]), P-2, P)) % P
        for i in range(n):
            if i != r and M[i,c]: M[i] = (M[i] - M[i,c]*M[r]) % P
        r += 1
        if r == n: break
    return r

# use ALL 99 Gargs as candidates (they are the allowed small list); class vectors
gamA = {g: intvec(cls235(col_Li3(sp.Rational(g), SS235))) for g in Gargs}
tau_i = {e: intvec(v) for e,v in tau.items()}

sol_per_edge = {}
for e in (1,2,3):
    te = tau_i[e]
    # k=1
    one = None
    for g, gv in gamA.items():
        if grank_rows([gv, te]) == 1:
            one = g; break
    if one is not None:
        sol_per_edge[e] = [one]; print(f"edge {e}: 1 column: {one}"); continue
    # k=2
    keys = list(gamA.keys())
    two = None
    for i in range(len(keys)):
        gi = gamA[keys[i]]
        for j in range(i+1, len(keys)):
            gj = gamA[keys[j]]
            if grank_rows([gi, gj]) == 2 and grank_rows([gi, gj, te]) == 2:
                two = (keys[i], keys[j]); break
        if two: break
    sol_per_edge[e] = list(two) if two else None
    print(f"edge {e}: 2 columns: {two}   ({time.time()-t0:.1f}s)")
print(sol_per_edge)
