import sys, os
sys.dont_write_bytecode = True
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
from towlib import *
from genintegrand import S_poly
from fractions import Fraction as Fr

tstar = QSI(Fr(1,2),0,Fr(1,2),0)
def S_val(d, m, l):
    P = S_poly(d, m)
    out = ZERO
    for n, c in enumerate(P):
        if n < l: continue
        ff = 1
        for a in range(l): ff *= (n-a)
        out = out + QSI(Fr(c*ff)) * (tstar**(n-l))
    return out

# complex number as (re, im) Fraction pair; Gaussian elimination over Q(i)
def solve_QI(A, b):
    # A: list of rows of (re,im) pairs; b: list of pairs. Returns solution list of pairs.
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = None
        for r in range(col, n):
            if M[r][col] != (Fr(0), Fr(0)):
                piv = r; break
        assert piv is not None, 'singular'
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        inv = (pv[0]/(pv[0]**2+pv[1]**2), -pv[1]/(pv[0]**2+pv[1]**2))
        M[col] = [((c[0]*inv[0]-c[1]*inv[1]), (c[0]*inv[1]+c[1]*inv[0])) for c in M[col]]
        for r in range(n):
            if r != col and M[r][col] != (Fr(0),Fr(0)):
                f = M[r][col]
                M[r] = [ (c[0]-f[0]*m[0]+f[1]*m[1], c[1]-f[0]*m[1]-f[1]*m[0]) for c,m in zip(M[r], M[col]) ]
    return [M[i][n] for i in range(n)]

def fitR(l, dmax=12):
    deg = 2*l
    mons = [(a,b) for a in range(0,deg+1) for b in range(0,deg+1-a)]
    n = len(mons)
    pts = []
    for d in range(l+1, dmax+1):
        for m in range(0, d+1):
            pts.append((d,m))
    # pick a square subsystem greedily to keep numbers small, then verify on all
    chosen = []
    rowsC = []
    def row_of(pt):
        d,m = pt
        return [(Fr(d**a*m**b), Fr(0)) for (a,b) in mons]
    import itertools
    # simple approach: use normal equations with Fraction arithmetic
    G = [[(Fr(0),Fr(0)) for _ in range(n)] for _ in range(n)]
    gb = [(Fr(0),Fr(0)) for _ in range(n)]
    for pt in pts:
        d,m = pt
        R = S_val(d,m,l)/S_val(d,m,0)
        prod = 1
        for j in range(1,l+1): prod *= (2*d-2*j-1)
        val = (R.p0*prod, R.p2*prod)
        row = row_of(pt)
        for a in range(n):
            for bq in range(n):
                G[a][bq] = (G[a][bq][0]+row[a][0]*row[bq][0], G[a][bq][1])
            gb[a] = (gb[a][0]+row[a][0]*val[0], gb[a][1]+row[a][0]*val[1])
    sol = solve_QI(G, gb)
    # verify
    maxres = Fr(0)
    for pt in pts:
        d,m = pt
        R = S_val(d,m,l)/S_val(d,m,0)
        prod = 1
        for j in range(1,l+1): prod *= (2*d-2*j-1)
        val = (R.p0*prod, R.p2*prod)
        row = row_of(pt)
        acc = (Fr(0),Fr(0))
        for a in range(n):
            acc = (acc[0]+row[a][0]*sol[a][0], acc[1]+row[a][0]*sol[a][1])
        maxres = max(maxres, abs(acc[0]-val[0]), abs(acc[1]-val[1]))
    return mons, sol, maxres, len(pts)

import json
out = {}
for l in range(1,5):
    mons, sol, res, npts = fitR(l)
    print(f'=== l={l}, {npts} pts, max residual {res}')
    terms = []
    for (a,b), (re,im) in zip(mons, sol):
        if re or im:
            terms.append(f'({"+" if im>=0 else ""}{im}i{"+" if re>=0 and im==0 else ""}{re if re or not im else ""})*d^{a}m^{b}' if (a or b) else f'({re}{"+" if im>=0 else ""}{im}i)')
    # prettier: print as dict
    out[l] = [(list(ab), str(re), str(im)) for ab,(re,im) in zip(mons,sol) if re or im]
    for ab,(re,im) in zip(mons,sol):
        if re or im:
            print(f'   d^{ab[0]} m^{ab[1]}: {re} {"+" if im>=0 else ""}{im}*i')
json.dump(out, open('Sderiv_fits.json','w'), indent=1)
