script = r'''"""M4e: value-relations among the symbol-zero columns [kv_indep, L, u_z] (no target).

Determines the solution-coset dimension and whether L(chi8,5) is entangled in
relations (making its coefficient representative-dependent).
"""
import json, math, signal, sys, time
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
def lll_rows(vals, digits):
    n=len(vals); C=mp.mpf(10)**digits
    rows=[[int(i==j) for j in range(n)]+[int(mp.nint(C*v))] for i,v in enumerate(vals)]
    signal.alarm(2700); red=flint.fmpz_mat(rows).lll(gram='exact'); signal.alarm(0)
    out=[]
    for i in range(n):
        vec=[int(red[i,j]) for j in range(n)]
        out.append(dict(vec=vec,maxcoef=max(map(abs,vec)),
                        residual=mp.fsum(c*v for c,v in zip(vec,vals))))
    return out
def solve_for_target(basis, cols, t):
    words = NS.WORDS
    tp=[F(str(t.get(w,ZERO).p)) for w in words]; tq=[F(str(t.get(w,ZERO).q)) for w in words]
    rows=[[str(c) for c in r]+[str(tp[i]),str(tq[i])] for i,r in enumerate(zip(*cols))]
    mat,rank=flint.fmpq_mat(rows).rref()
    n=len(basis)
    pivots=[next(j for j in range(mat.ncols()) if mat[i,j]) for i in range(rank)]
    xp=[F(0)]*n; xq=[F(0)]*n
    for i,p in enumerate(pivots): xp[p]=F(str(mat[i,n])); xq[p]=F(str(mat[i,n+1]))
    return xp,xq
def main(digits=1200):
    ctx.prec = math.ceil((digits+140)*math.log2(10)); mp.mp.dps = digits+120
    S=list(NS.S0); sol=NS.solve(S,'vr'); basis=sol['basis']; cols=NS.symbol_columns(basis)
    ks=json.loads((HERE/'nice_ksat_s0.json').read_text()); K=ks['K']
    recipes=[(b[0],b[2]) for b in basis]
    vals=[mp.mpf(v.mid().str(digits+110,radius=False)) for v in [evaluate(r,digits+60) for _,r in recipes]]
    r2=mp.sqrt(2)
    Lchi=mp.mpf(json.loads((HERE/'Lchi8_5.json').read_text())['value'])
    kv=[mp.fsum(F(c)*v for c,v in zip(krow,vals)) for krow in K]
    dropk={17,27,30}; kindep=[j for j in range(35) if j not in dropk]
    zlist=[NS.EXTRA[i] for i in [0,2,3,4,5]]
    uzs=[]
    for z in zlist:
        xpz,xqz=solve_for_target(basis,cols,W.li_symbol(5,z))
        pv=mp.fsum((mp.mpf(x.numerator)/x.denominator+r2*mp.mpf(q.numerator)/q.denominator)*v for x,q,v in zip(xpz,xqz,vals))
        zval=mp.mpf(evaluate(('li',5,z),digits+60).mid().str(digits+110,radius=False))
        uzs.append(zval-pv)
    cn=[f'kv{j}' for j in kindep]+[f's2kv{j}' for j in kindep]+['L','s2L']+[f'u{k}' for k in range(5)]+[f's2u{k}' for k in range(5)]
    colv=[kv[j] for j in kindep]+[r2*kv[j] for j in kindep]+[Lchi,r2*Lchi]+uzs+[r2*u for u in uzs]
    rows=lll_rows(colv,digits)
    tol=mp.mpf(10)**(-digits-20)
    deps=[r for r in rows if abs(r['residual'])<tol]
    deps.sort(key=lambda r:r['maxcoef'])
    log(f'{len(colv)} columns; value-relations found: {len(deps)}')
    for r in deps:
        sup=[(cn[j],c) for j,c in enumerate(r['vec']) if c]
        log(f"  maxcoef={r['maxcoef']} resid={mp.nstr(r['residual'],4)}")
        log(f"    {sup}")
if __name__=='__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 1200)
'''
M4E = '<project>/research/08_theory/M4e_w5_lchi'
open(M4E+'/lchi_valrel.py','w').write(script)
import subprocess, time
PY = '<project>/.venv/bin/python'
t0=time.monotonic()
r = subprocess.run(['nice', PY, 'lchi_valrel.py', '1200'], cwd=M4E, capture_output=True, text=True, timeout=2700)
print(r.stdout); print("STDERR:", r.stderr[-800:]); print(f"elapsed {time.monotonic()-t0:.0f}s")
open(M4E+'/progress.log','a').write("step_56: lchi_valrel.py 1200 (value-relation structure)\n")
