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
gamA = {g: cls235(col_Li3(g, SS235)) for g in poolA}

# integer-scaled class vectors
def intvec(v):
    den = sp.ilcm(*[x.q for x in v])
    return tuple(int(sp.Rational(x)*den) for x in v)
tau_i = {e: intvec(v) for e,v in tau.items()}
gam_i = {g: intvec(v) for g,v in gamA.items()}
Tv = np.array([tau_i[1], tau_i[2], tau_i[3]], dtype=object)  # exact ints

P = 2**31-1
def modmat(vecs):
    return np.array([[x % P for x in v] for v in vecs], dtype=np.int64)
def grank(M):
    M = M.copy() % P
    n, m = M.shape
    r = 0
    for c in range(m):
        # find pivot
        piv = -1
        for i in range(r, n):
            if M[i,c] % P: piv = i; break
        if piv < 0: continue
        if piv != r: M[[r,piv]] = M[[piv,r]]
        inv = pow(int(M[r,c]), P-2, P)
        M[r] = (M[r]*inv) % P
        for i in range(n):
            if i != r and M[i,c]:
                M[i] = (M[i] - M[i,c]*M[r]) % P
        r += 1
        if r == n: break
    return r

Tmod = modmat([tau_i[1], tau_i[2], tau_i[3]])  # 3x8
rkT = grank(Tmod); print("rank T mod p =", rkT)
keys = list(gam_i.keys())
Gmod = {g: modmat([gam_i[g]])[0] for g in keys}  # 8-vectors

# search quadruples: enumerate triples with dimcap==2, rank==3, then 4th
def dimcap_mod(sel):
    if not sel: return 0
    S = modmat([gam_i[g] for g in sel])
    rS = grank(S)
    ST = np.vstack([S, Tmod])
    return rkT + rS - grank(ST), rS
trip_count = 0
found = None
t1=time.time()
for trip in itertools.combinations(range(len(keys)), 3):
    sel = [keys[i] for i in trip]
    d, rS = dimcap_mod(sel)
    if d == 2 and rS == 3:
        # try all fourth candidates
        S = modmat([gam_i[g] for g in sel])
        ST = np.vstack([S, Tmod])
        base_r = grank(ST)
        for g4 in keys:
            if g4 in sel: continue
            S4 = np.vstack([S, Gmod[g4][None,:]])
            r4 = grank(S4)
            if r4 != 4: continue
            r4T = grank(np.vstack([S4, Tmod]))
            if r4T == 4:
                found = sel + [g4]; break
        if found: break
    trip_count += 1
print("searched", trip_count, "triples; found quadruple:", found, f"({time.time()-t1:.1f}s)")
if found:
    # verify over Q
    S = sp.Matrix.hstack(*[gamA[g] for g in found])
    ST = sp.Matrix.hstack(S, tau[1], tau[2], tau[3])
    print("over Q: rank S =", S.rank(), "rank [S|T] =", ST.rank(), " -> T subset span:", S.rank()==ST.rank()==4)
    pickle.dump({'chosen123': found}, open("short123.pkl","wb"))
