WORK='<project>/research/08_theory/M4f_w5_extras'
src = '''"""M4f step 1: pool of exceptional S-units x in (-1,1)\\{0} of Z[sqrt2,1/2],
x and 1-x both +- (sqrt2)^m (1+sqrt2)^n, |m|<=12, |n|<=8.
For each: does li_symbol(5,x) lie in the Q-span of the S0 atom symbols (nice_search.S0)?
If yes, record the exact particular solution p_x over the 65 S0 atoms (needed for u_x).
Saves pool.json, ordered by height h = max(|mx|+|nx|, |m1x|+|n1x|).
"""
import json, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
import w5_symbol as W
import nice_search as NS
from qsqrt2 import sunit_log

ZERO, ONE, Q, R, A, B = W.ZERO, W.ONE, W.Q, W.R, W.A, W.B
def log(s): print(s, flush=True)

def solve_for_target(basis, cols, t):
    words = NS.WORDS
    tp=[F(str(t.get(w,ZERO).p)) for w in words]; tq=[F(str(t.get(w,ZERO).q)) for w in words]
    rows=[[str(c) for c in r]+[str(tp[i]),str(tq[i])] for i,r in enumerate(zip(*cols))]
    mat,rank=flint.fmpq_mat(rows).rref()
    n=len(basis)
    pivots=[next(j for j in range(mat.ncols()) if mat[i,j]) for i in range(rank)]
    if not all(j<n for j in pivots): return None
    xp=[F(0)]*n; xq=[F(0)]*n
    for i,p in enumerate(pivots): xp[p]=F(str(mat[i,n])); xq[p]=F(str(mat[i,n+1]))
    return xp,xq

def main():
    # 1. enumerate
    raw = {}
    for eps in [1,-1]:
        for m in range(-12,13):
            Rm = R**m
            for n in range(-8,9):
                x = Rm*(A**n)
                if eps<0: x = -x
                raw[x] = (eps,m,n)
    log(f'candidate S-units: {len(raw)}')
    pool = {}
    for x,(eps,m,n) in raw.items():
        if not (-ONE < x < ONE) or x == ZERO: continue
        try:
            jx,kx = sunit_log(x)     # |x| = R^jx A^kx
            jy,ky = sunit_log(ONE-x) # |1-x| = R^jy A^ky
        except ValueError:
            continue
        if max(abs(jx),abs(jy))>12 or max(abs(kx),abs(ky))>8: continue
        sx = x.sign(); sy = (ONE-x).sign()
        pool[x] = dict(mx=sx*jx, nx=sx*kx*0+kx, m1x=jy, n1x=ky,
                       signx=sx, sign1x=sy,
                       h=max(abs(jx)+abs(kx), abs(jy)+abs(ky)))
    log(f'exceptional units in (-1,1): {len(pool)}')
    # 2. symbol membership in S0 span
    sol = NS.solve(list(NS.S0), 'pool')
    basis = sol['basis']; cols = NS.symbol_columns(basis)
    keep = []
    for x in sorted(pool, key=lambda z:(pool[z]['h'], abs(float(z)))):
        t = W.li_symbol(5, x)
        out = solve_for_target(basis, cols, t)
        d = pool[x]
        line = (f'x={str(x):24s} |x|=R^{d["mx"]}A^{d["nx"]} |1-x|=R^{d["m1x"]}A^{d["n1x"]} '
                f'h={d["h"]:2d} sign={d["signx"]:+d}')
        if out is None:
            log(line + '  symbol NOT in span')
            continue
        xp,xq = out
        dn = 1
        for c in xp+xq: dn = dn*c.denominator//__import__("math").gcd(dn,c.denominator)
        mx = max([abs(c.numerator) for c in xp+xq]+[0])
        log(line + f'  IN SPAN; particular max|num|={mx} lcm-den={dn}')
        keep.append(dict(p=str(x.p), q=str(x.q), h=d['h'],
                         mx=d['mx'], nx=d['nx'], m1x=d['m1x'], n1x=d['n1x'], signx=d['signx'],
                         xp=[str(c) for c in xp], xq=[str(c) for c in xq],
                         pmax=mx, pden=dn))
    json.dump(dict(names=sol['names'], pool=keep), open(HERE/'pool.json','w'), indent=1)
    log(f'pool (symbol in span): {len(keep)}; saved pool.json')

if __name__=='__main__':
    main()
'''
open(f'{WORK}/pool.py','w').write(src)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','pool.py'],
                   cwd=WORK, capture_output=True, text=True, timeout=1500)
print(r.stdout)
print('STDERR', r.stderr[-2000:])
