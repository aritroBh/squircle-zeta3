import sys, time, itertools, pickle
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
def intvec(v):
    den = sp.ilcm(*[x.q for x in v])
    return np.array([int(sp.Rational(x)*den) for x in v], dtype=object)
gamA = {g: intvec(cls235(col_Li3(sp.Rational(g), SS235))) for g in Gargs}
gamA = {g: v for g, v in gamA.items() if any(v)}
tau_i = {e: intvec(v) for e,v in tau.items()}
keys = list(gamA.keys())
print("candidates:", len(keys))

P = 2**31-1
G = np.array([[int(x) % P for x in gamA[g]] for g in keys], dtype=np.int64)  # 96 x 8
def rref_rows(M):
    M = M.copy() % P; n, m = M.shape; r = 0; pivs = []
    for c in range(m):
        piv = -1
        for i in range(r, n):
            if M[i,c]: piv = i; break
        if piv < 0: continue
        if piv != r: M[[r,piv]] = M[[piv,r]]
        M[r] = (M[r]*pow(int(M[r,c]), P-2, P)) % P
        for i in range(n):
            if i != r and M[i,c]: M[i] = (M[i] - M[i,c]*M[r]) % P
        pivs.append(c); r += 1
        if r == n: break
    return M[:r], pivs

def in_span_matrix(basis_rref, pivs, C):
    """basis_rref: r x m RREF rows with unit pivots at pivs; C: k x m. Returns boolean array: each row of C in span?"""
    if len(basis_rref) == 0: return np.zeros(C.shape[0], dtype=bool)
    R = C[:, pivs] % P          # k x r
    # candidate in span iff C - R*basis == 0
    resid = (C - R @ basis_rref) % P
    return (resid == 0).all(axis=1)

sols = {}
for e in (1,2,3):
    te = np.array([int(x) % P for x in tau_i[e]], dtype=np.int64)
    found_e = []
    t1 = time.time()
    for trip in itertools.combinations(range(len(keys)), 3):
        W3m = G[list(trip)]
        bW, pW = rref_rows(W3m)
        if len(bW) < 3: continue
        if in_span_matrix(bW, pW, te[None,:])[0]: continue  # k<=3 (none expected)
        bW4, pW4 = rref_rows(np.vstack([W3m, te]))
        hits = np.nonzero(in_span_matrix(bW4, pW4, G))[0]
        for h in hits:
            if h in trip: continue
            if not in_span_matrix(bW, pW, G[h][None,:])[0]:
                found_e.append((keys[h],) + tuple(keys[i] for i in trip))
        if len(found_e) >= 30: break
    print(f"edge {e}: found {len(found_e)} 4-column solutions ({time.time()-t1:.1f}s), first: {found_e[0] if found_e else None}")
    sols[e] = found_e
pickle.dump(sols, open("sols4_modp.pkl","wb"))
print(f"total ({time.time()-t0:.1f}s)")
