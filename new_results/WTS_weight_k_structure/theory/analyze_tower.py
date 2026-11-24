"""TOWER-LEVEL pattern: for every key (i,j) of degree n and every pole order r (2<=r<=2n-2) at v=+i,
the coefficient h_{ij}^{(k)}(r) in Q(sqrt2,i) as a function of the level k >= n.
Report the minimal polynomial degree in k (exact finite differences over k>=n+1, and also whether k=n fits).
Depth delta = (2n-2) - r.  Prints a table: for each (n, delta) the max degree over keys/orders,
and the first k0 from which the polynomial law holds."""
import sys, os, pickle, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
from fractions import Fraction as Fr
T = {}
for fn in sorted(glob.glob(os.path.join(HERE, "tower_*.pkl"))):
    T[int(os.path.basename(fn)[6:-4])] = pickle.load(open(fn, 'rb'))
ks = sorted(T); K = ks[-1]
print("levels:", ks)
def fd(seq, m):
    s = list(seq)
    for _ in range(m): s = [s[i+1]-s[i] for i in range(len(s)-1)]
    return s
def getc(k, key, r):
    return T[k].get(key, {}).get(r, (Fr(0),)*4)
def mindeg(seq):
    for m in range(len(seq)):
        if all(x == 0 for x in fd(seq, m)):
            return m-1
    return None
res = {}
keys = sorted({kk for k in ks for kk in T[k]}, key=lambda kk: (sum(kk), -kk[0]))
for key in keys:
    n = sum(key)
    if n < 2: continue
    for r in range(2, 2*n-1):
        delta = 2*n-2-r
        best = -1; start_ok = None
        for comp in range(4):
            seq = [getc(k, key, r)[comp] for k in range(n+1, K+1)]
            if len(seq) < 3:  # not testable
                best = None; break
            md = mindeg(seq)
            if md is None: best = None; break
            best = max(best, md)
        # does k=n lie on the polynomial?
        onpoly = None
        if best is not None and best >= 0 and (K - n) >= best + 3:
            onpoly = True
            for comp in range(4):
                seq = [getc(k, key, r)[comp] for k in range(n+1, K+1)]
                pts = list(zip(range(n+1, n+1+best+1), seq[:best+1]))
                val = Fr(0)
                for a_, (ka, ya) in enumerate(pts):
                    term = ya
                    for b_, (kb, yb) in enumerate(pts):
                        if a_ != b_: term *= Fr(n-kb, ka-kb)
                    val += term
                if val != getc(n, key, r)[comp]: onpoly = False
        res[(key, r)] = (n, delta, best, onpoly, len(range(n+1, K+1)))
# table by (n, delta)
import collections
tab = collections.defaultdict(list)
for (key, r), (n, delta, best, onpoly, npts) in res.items():
    tab[(n, delta)].append((key, r, best, onpoly, npts))
for (n, delta) in sorted(tab):
    L = tab[(n, delta)]
    degs = sorted({x[2] for x in L}, key=lambda z: (z is None, z))
    onp = sorted({x[3] for x in L}, key=lambda z: str(z))
    print(f"n={n} depth delta={delta} (order r={2*n-2-delta}): poly degree in k over keys = {degs};  k=n on polynomial: {onp}; npts={L[0][4]}")
