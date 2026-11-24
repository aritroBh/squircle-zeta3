"""M4e: minimal extras for target = S0-kernel + L + sum a_z u_z.

S0 particular (nice_symbol.json) + saturated S0 kernel (nice_ksat_s0.json).
u_z = Li5(z) - val(S0 particular for symbol Li5(z)), symbol-zero, for the 5
extras that are NOT already in span(V0) (skip b^4).
Dependent S0 kernel columns dropped upfront (kv17, kv27, kv30 + sqrt2).
Reports which u_z are used and the relation.
"""
import json, math, pickle, signal, sys, time
from fractions import Fraction as F
from math import gcd, lcm
from pathlib import Path
import flint
from flint import arb, ctx
import mpmath as mp
import w5_symbol as W
from w5_search import evaluate, mpl, numeric
import nice_search as NS

HERE = Path(__file__).resolve().parent
ZERO, ONE, Q, R, A, B = W.ZERO, W.ONE, W.Q, W.R, W.A, W.B
def log(s): print(s, flush=True)
def qnum(c): return mp.mpf(c.p.numerator)/c.p.denominator + mp.sqrt(2)*mp.mpf(c.q.numerator)/c.q.denominator
def lll_rows(vals, digits):
    n=len(vals); C=mp.mpf(10)**digits
    rows=[[int(i==j) for j in range(n)]+[int(mp.nint(C*v))] for i,v in enumerate(vals)]
    signal.alarm(2700); t0=time.monotonic()
    red=flint.fmpz_mat(rows).lll(gram='exact'); signal.alarm(0)
    out=[]
    for i in range(n):
        vec=[int(red[i,j]) for j in range(n)]
        residual=mp.fsum(c*v for c,v in zip(vec,vals))
        out.append(dict(vec=vec,maxcoef=max(map(abs,vec)),residual=residual))
    log(f'  LLL {n}x{n}@{digits}: {time.monotonic()-t0:.1f}s')
    return out
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

def main(digits=1200):
    ctx.prec = math.ceil((digits+140)*math.log2(10)); mp.mp.dps = digits+120
    S = list(NS.S0)
    sol = NS.solve(S, 'min')
    names = sol['names']; basis = sol['basis']; xp0, xq0 = sol['xp'], sol['xq']
    cols = NS.symbol_columns(basis)
    ks = json.loads((HERE/'nice_ksat_s0.json').read_text()); K = ks['K']
    recipes = [(b[0], b[2]) for b in basis]
    values = [evaluate(r, digits+60) for _, r in recipes]
    vals = [mp.mpf(v.mid().str(digits+110, radius=False)) for v in values]
    r2 = mp.sqrt(2)
    Lchi = mp.mpf(json.loads((HERE/'Lchi8_5.json').read_text())['value'])
    cert = json.loads((HERE/'cobracket.json').read_text())
    fixedv = mp.mpf(0)
    for term in cert['correction']:
        x,y,c = (Q(*map(F,term[k])) for k in ['x','y','coefficient'])
        v = evaluate(('mpl',3,2,x,y), digits+60)
        fixedv += qnum(c)*mp.mpf(v.mid().str(digits+110, radius=False))
    saved = json.loads((HERE/'W5_5400.json').read_text())
    pval = mp.fsum((mp.mpf(x.numerator)/x.denominator + r2*mp.mpf(q.numerator)/q.denominator)*v
                   for x,q,v in zip(xp0,xq0,vals))
    target = mp.mpf(saved['value']) - fixedv - pval
    kv = [mp.fsum(F(c)*v for c,v in zip(krow,vals)) for krow in K]
    # independent S0 kernel columns: drop kv17, kv27, kv30
    dropk = {17, 27, 30}
    kindep = [j for j in range(35) if j not in dropk]
    # u_z for 5 extras (skip b^4 = EXTRA[1])
    zlist = [NS.EXTRA[i] for i in [0,2,3,4,5]]
    uzs = []
    for z in zlist:
        t = W.li_symbol(5, z)
        out = solve_for_target(basis, cols, t)
        xpz, xqz = out
        pv = mp.fsum((mp.mpf(x.numerator)/x.denominator + r2*mp.mpf(q.numerator)/q.denominator)*v
                     for x,q,v in zip(xpz,xqz,vals))
        zval = mp.mpf(evaluate(('li',5,z), digits+60).mid().str(digits+110, radius=False))
        uzs.append(zval - pv)
    cn = ['target'] + [f'kv{j}' for j in kindep] + [f's2kv{j}' for j in kindep] + \
         ['L','s2L'] + [f'u{z}' for z in zlist] + [f's2u{z}' for z in zlist]
    colv = [target] + [kv[j] for j in kindep] + [r2*kv[j] for j in kindep] + \
           [Lchi, r2*Lchi] + uzs + [r2*u for u in uzs]
    tol = mp.mpf(10)**(-digits-20)
    rows = lll_rows(colv, digits)
    gen = [r for r in rows if r['vec'][0]!=0 and abs(r['residual'])<tol]
    gen.sort(key=lambda r: r['maxcoef'])
    if not gen:
        log('NO relation'); return
    r = gen[0]; vec=r['vec']; c0=vec[0]
    log(f'c0={c0} maxcoef={r["maxcoef"]} resid={mp.nstr(r["residual"],6)}')
    nu = len(zlist)
    for k,z in enumerate(zlist):
        a = vec[1+len(kindep)*2+2+k]; b = vec[1+len(kindep)*2+2+nu+k]
        if a or b:
            log(f'  uses u({z}): rat {F(-a,c0)}, s2 {F(-b,c0)}')
    log(f'  L: {F(-vec[1+len(kindep)*2],c0)} + {F(-vec[1+len(kindep)*2+1],c0)} sqrt2')
    json.dump(dict(digits=digits,c0=c0,cn=cn,vec=vec,resid=mp.nstr(r['residual'],8)),
              open(HERE/f'lchi_min_{digits}.json','w'), indent=1)

if __name__=='__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 1200)
