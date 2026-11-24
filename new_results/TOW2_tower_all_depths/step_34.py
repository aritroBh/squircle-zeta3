import sys, os, json
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

# refit l=4 deterministically (same as before but print correctly)
def solve_QI(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = None
        for r in range(col, n):
            if M[r][col] != (Fr(0), Fr(0)): piv = r; break
        if piv is None: return None
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        inv = (pv[0]/(pv[0]**2+pv[1]**2), -pv[1]/(pv[0]**2+pv[1]**2))
        M[col] = [((c[0]*inv[0]-c[1]*inv[1]), (c[0]*inv[1]+c[1]*inv[0])) for c in M[col]]
        for r in range(n):
            if r != col and M[r][col] != (Fr(0),Fr(0)):
                f = M[r][col]
                M[r] = [ (c[0]-f[0]*m[0]+f[1]*m[1], c[1]-f[0]*m[1]-f[1]*m[0]) for c,m in zip(M[r], M[col]) ]
    return [M[i][n] for i in range(n)]

deg = 8; l = 4
mons = [(a,b) for a in range(0,deg+1) for b in range(0,deg+1-a)]
n = len(mons)
pts = [(d,m) for d in range(l+1, 15) for m in range(0, d+1)]
def row_of(pt):
    d,m = pt
    return [(Fr(d**a*m**b), Fr(0)) for (a,b) in mons]
def val_of(pt):
    d,m = pt
    R = S_val(d,m,l)/S_val(d,m,0)
    prod = 1
    for j in range(1,l+1): prod *= (2*d-2*j-1)
    return (R.p0*prod, R.p2*prod)
sub = pts[::3][:n]
if len(sub) < n: sub = pts[:n]
sol = solve_QI([row_of(p) for p in sub], [val_of(p) for p in sub])
assert sol is not None
ok = True
for p in pts:
    row = row_of(p); val = val_of(p)
    acc = (Fr(0),Fr(0))
    for a in range(n):
        acc = (acc[0]+row[a][0]*sol[a][0], acc[1]+row[a][0]*sol[a][1])
    if acc != val: ok = False; print('MISMATCH at', p); break
print('N4 fit verifies on all', len(pts), 'points:', ok)
terms = [(ab, re, im) for ab,(re,im) in zip(mons,sol) if re or im]
for ab,(re,im) in terms:
    print(f'   d^{ab[0]} m^{ab[1]}: {re} {"+" if im>=0 else ""}{im}*i')
json.dump([(list(ab),str(re),str(im)) for ab,(re,im) in terms], open('N4_fit.json','w'), indent=1)
