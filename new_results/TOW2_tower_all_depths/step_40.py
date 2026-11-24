import sys, os, json, random
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
from towlib import *
from fractions import Fraction as Fr
from math import comb

def e_from_data(key, delta, s):
    i,j = key; dd = i+j
    base = (2*i+1) if i==j else dd
    return sum( u(i,j,base+t, delta).scale(Fr((-1)**(s-t)*comb(s,t))) for t in range(0,s+1) )

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

random.seed(7)
def fit_e(delta, s, deg, den_pow2, den_pi, dmax=11):
    """fit e(delta,s)/e00 = N(d,m)/[2^den_pow2 * prod_{j=1..den_pi}(2d-2j-1)]"""
    mons = [(a,b) for a in range(0,deg+1) for b in range(0,deg+1-a)]
    n = len(mons)
    pts = []
    for dd in range(max(2, delta//2+2), dmax+1):
        for key in keys_of_degree(dd):
            i,j = key
            if delta > 2*dd-4: continue
            base = (2*i+1) if i==j else dd
            if base+s > 14: continue
            pts.append(key)
    def val_of(key):
        i,j = key; dd = i+j
        e = e_from_data(key, delta, s)
        e00 = e_from_data(key, 0, 0)
        R = e/e00
        den = Fr(2**den_pow2)
        for jj in range(1, den_pi+1): den *= (2*dd-2*jj-1)
        return (R.p0*den, R.p2*den)
    def row_of(key):
        dd = key[0]+key[1]; mm = key[0]
        return [(Fr(dd**a*mm**b), Fr(0)) for (a,b) in mons]
    for attempt in range(400):
        sub = random.sample(pts, n)
        sol = solve_QI([row_of(p) for p in sub], [val_of(p) for p in sub])
        if sol is None: continue
        ok = all(sum_c := True for _ in [0])
        ok = True
        for p in pts:
            row = row_of(p); val = val_of(p)
            acc = (Fr(0),Fr(0))
            for a2 in range(n):
                acc = (acc[0]+row[a2][0]*sol[a2][0], acc[1]+row[a2][0]*sol[a2][1])
            if acc != val: ok = False; break
        if ok:
            return mons, sol, len(pts)
    return None, None, len(pts)

# fit all (delta, s) with delta <= 6
table = {}
for delta in range(0, 7):
    for s in range(0, delta//2+1):
        if (delta, s) in [(0,0),(1,0),(2,0),(2,1),(3,0),(3,1),(4,0),(4,2),(5,2),(6,3)]:
            continue  # known or rho
        got = None
        for deg in (2, 3, 4, 5, 6):
            for den_pow2 in range(0, 7):
                for den_pi in range(0, delta+2):
                    mons, sol, npts = fit_e(delta, s, deg, den_pow2, den_pi)
                    if sol is not None:
                        got = (deg, den_pow2, den_pi, mons, sol, npts)
                        break
                if got: break
            if got: break
        if got:
            deg, den_pow2, den_pi, mons, sol, npts = got
            table[(delta,s)] = got
            nz = [(ab,re,im) for ab,(re,im) in zip(mons,sol) if re or im]
            print(f'e({delta},{s}): deg={deg} den=2^{den_pow2}*Pi_{den_pi}, {npts} pts, {len(nz)} terms')
        else:
            print(f'e({delta},{s}): NO FIT FOUND ({npts} pts)')
import pickle
pickle.dump(table, open('e_fits.pkl','wb'))
