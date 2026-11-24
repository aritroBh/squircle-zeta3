import sys, os, json
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
from towlib import *
from fractions import Fraction as Fr
from math import comb, gcd

def e_from_data(key, delta, s):
    i,j = key; dd = i+j
    base = (2*i+1) if i==j else dd
    return sum( u(i,j,base+t, delta).scale(Fr((-1)**(s-t)*comb(s,t))) for t in range(0,s+1) )

def Rvals(delta, s, dmax=11):
    out = []
    for dd in range(max(2, delta//2+2), dmax+1):
        for key in keys_of_degree(dd):
            i,j = key
            if delta > 2*dd-4: continue
            base = (2*i+1) if i==j else dd
            if base+s > 14: continue
            e = e_from_data(key, delta, s); e00 = e_from_data(key, 0, 0)
            R = e/e00
            out.append((key, (R.p0, R.p2)))
    return out

def den_structure(vals, maxb=10):
    # for each d: lcm of denominators; then find minimal (a,b) with lcm_d | 2^a Pi_b(d) for all d
    from collections import defaultdict
    byd = defaultdict(list)
    for key, (re, im) in vals:
        byd[key[0]+key[1]].extend([re.denominator, im.denominator])
    lcms = {d: 1 for d in byd}
    for d, ds in byd.items():
        L = 1
        for x in ds: L = L*x//gcd(L,x)
        lcms[d] = L
    def Pi(d, b):
        out = 1
        for j in range(1, b+1): out *= (2*d-2*j-1)
        return out
    for b in range(0, maxb+1):
        for a in range(0, 12):
            if all((2**a * Pi(d, b)) % L == 0 for d, L in lcms.items()):
                return a, b, lcms
    return None, None, lcms

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

def fit_e(delta, s, dmax=11):
    vals = Rvals(delta, s, dmax)
    a, b, lcms = den_structure(vals)
    if a is None: return None
    def den(dd):
        out = Fr(2**a)
        for j in range(1, b+1): out *= (2*dd-2*j-1)
        return out
    for deg in range(1, 11):
        mons = [(x,y) for x in range(0,deg+1) for y in range(0,deg+1-x)]
        n = len(mons)
        rows = []; rhs = []
        for key, (re, im) in vals:
            dd = key[0]+key[1]; mm = key[0]
            rows.append([(Fr(dd**x*mm**y), Fr(0)) for (x,y) in mons])
            rhs.append((re*den(dd), im*den(dd)))
        # greedy independent selection
        idx = list(range(len(rows)))
        sel = []
        cur = []
        # try in order of small d first
        for t in idx:
            trial = cur + [t]
            A = [rows[tt] for tt in trial]
            B = [rhs[tt] for tt in trial]
            if len(trial) < n:
                # rank check via attempting elimination on extended basis
                pass
            sel.append(t)
            if len(sel) == n: break
        sol = solve_QI([rows[t] for t in sel], [rhs[t] for t in sel])
        if sol is None: continue
        ok = True
        for t in idx:
            acc = (Fr(0),Fr(0))
            for q in range(n):
                acc = (acc[0]+rows[t][q][0]*sol[q][0], acc[1]+rows[t][q][0]*sol[q][1])
            if acc != rhs[t]: ok = False; break
        if ok:
            return dict(delta=delta, s=s, a=a, b=b, deg=deg, mons=mons, sol=sol, npts=len(vals), lcms={str(k):v for k,v in lcms.items()})
    return None

results = {}
for (delta, s) in [(4,1),(5,0),(5,1),(6,0),(6,1),(6,2)]:
    r = fit_e(delta, s)
    if r is None:
        print(f'e({delta},{s}): NO FIT')
    else:
        results[(delta,s)] = r
        print(f'e({delta},{s}): den=2^{r["a"]}*Pi_{r["b"]}, deg={r["deg"]}, {r["npts"]} pts OK')
import pickle
pickle.dump(results, open('e_fits.pkl','wb'))
