import sys, time, itertools, pickle
sys.path.insert(0, ".")
from oct9d_common import *
import flint
import numpy as np
t0=time.time()
TW = TW3()
basis0 = [sp.Rational(g) for g in pickle.load(open("basis0.pkl","rb"))['basis0']]

def to_fmpq_mat(cols):
    """cols: list of sympy column matrices (rationals) -> flint fmpq_mat"""
    n = cols[0].rows
    M = flint.fmpq_mat(n, len(cols))
    for j, c in enumerate(cols):
        for i in range(n):
            x = sp.Rational(c[i])
            M[i, j] = flint.fmpq(x.p, x.q)
    return M

Gsym0 = {g: col_Li3(sp.Rational(g), SS2357) for g in sols7}
G20cols = [Gsym0[g] for g in basis0]

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

# mod-p rank machinery for greedy basis construction
P = 2**31 - 1
def intcol(c):
    den = sp.ilcm(*[sp.Rational(x).q for x in c])
    return np.array([int(sp.Rational(x)*den) % P for x in c], dtype=np.int64)
def grank_cols(cols):
    if not cols: return 0
    M = np.stack(cols, axis=1) % P
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

G20i = [intcol(c) for c in G20cols]
cur = list(G20i)
PB0 = []
for w in worder:
    for p in (2,3,5,7):
        c = col_prod(w, p, SS2357)
        ci = intcol(c)
        if grank_cols(cur + [ci]) > len(cur):
            PB0.append((w,p)); cur.append(ci)
        if len(cur) == 64: break
    if len(cur) == 64: break
print("PB0 count:", len(PB0), f"({time.time()-t0:.1f}s)")
print("PB0:", [(str(w),p) for w,p in PB0])
pickle.dump({'PB0': [(str(w),p) for w,p in PB0]}, open("PB0.pkl","wb"))

# exact solve with flint
cols = G20cols + [col_prod(w, p, SS2357) for w,p in PB0]
Mf = to_fmpq_mat(cols)
tf = flint.fmpq_mat(64, 1)
for i in range(64):
    x = sp.Rational(TW[0][i]); tf[i,0] = flint.fmpq(x.p, x.q)
print("flint det...")
d = Mf.det(); print("det =", d, f"({time.time()-t0:.1f}s)")
xf = Mf.solve(tf)
# verify
print("solve ok:", Mf*xf == tf, f"({time.time()-t0:.1f}s)")
sol0 = []
for i in range(64):
    q = xf[i,0]
    sol0.append(sp.Rational(sp.Integer(int(q.numerator)), sp.Integer(int(q.denominator))))
cl = sol0[:20]; cp = sol0[20:]
print("Li3 coeffs:", [str(c) for c in cl])
print("prod coeffs:", [str(c) for c in cp])
pickle.dump({'sol0': [str(c) for c in sol0]}, open("sol0_short.pkl","wb"))
