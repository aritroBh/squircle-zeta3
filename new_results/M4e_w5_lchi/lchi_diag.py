"""M4e diagnostic: the missing weight-5 direction.

Extras setup S0+EXTRA[:6] (71 atoms, kernel 41). Three LLL tests at `digits`:
 A: [target_e, kve, s2kve]            -> reproduce M4d's relation (vec0 ~ 1e8)
 B: [target_e, kve, s2kve, L, s2L]    -> small relation iff L(chi8,5) explains extras
 C: [L, kve, s2kve]                   -> is L in the extras-kernel value span?
Also prints the quotient structure.
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
def prim(vec):
    D = lcm(*[c.denominator for c in vec]); ints=[int(c*D) for c in vec]; g=0
    for i in ints: g=gcd(g,abs(i))
    return [i//g for i in ints]
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

def main(digits=800):
    ctx.prec = math.ceil((digits+140)*math.log2(10)); mp.mp.dps = digits+120
    S = list(NS.S0) + NS.EXTRA[:6]
    sol = NS.solve(S, 'diag')
    assert sol['ok']
    names = sol['names']; xp, xq = sol['xp'], sol['xq']
    K = [prim(v) for v in sol['kernel']]
    red = flint.fmpz_mat(K).lll(gram='exact')
    K = [[int(red[i,j]) for j in range(red.ncols())] for i in range(red.nrows())]
    log(f'atoms={len(names)} kernel={len(K)}')
    recipes = [(b[0], b[2]) for b in sol['basis']]
    t0=time.monotonic()
    values = [evaluate(r, digits+60) for _, r in recipes]
    vals = [mp.mpf(v.mid().str(digits+110, radius=False)) for v in values]
    log(f'atom values {time.monotonic()-t0:.1f}s')
    r2 = mp.sqrt(2)
    cert = json.loads((HERE/'cobracket.json').read_text())
    fixedv = mp.mpf(0)
    for term in cert['correction']:
        x,y,c = (Q(*map(F,term[k])) for k in ['x','y','coefficient'])
        v = evaluate(('mpl',3,2,x,y), digits+60)
        fixedv += qnum(c)*mp.mpf(v.mid().str(digits+110, radius=False))
    saved = json.loads((HERE/'W5_5400.json').read_text())
    pval = mp.fsum((mp.mpf(x.numerator)/x.denominator + r2*mp.mpf(q.numerator)/q.denominator)*v
                   for x,q,v in zip(xp,xq,vals))
    target = mp.mpf(saved['value']) - fixedv - pval
    kv = [mp.fsum(F(c)*v for c,v in zip(krow,vals)) for krow in K]
    Lchi = mp.mpf(json.loads((HERE/'Lchi8_5.json').read_text())['value'])
    tol = mp.mpf(10)**(-digits-20)
    def report(rows, label, vital0=True):
        gen = [r for r in rows if (r['vec'][0]!=0)==vital0 and abs(r['residual'])<tol]
        gen.sort(key=lambda r: r['maxcoef'])
        log(f'  [{label}] found {len(gen)}')
        for r in gen[:4]:
            log(f"    vec0={r['vec'][0]} maxcoef={r['maxcoef']} resid={mp.nstr(r['residual'],4)}")
        return gen
    colsA = [target]+kv+[r2*v for v in kv]
    log('TEST A: target in extras-kernel span (M4d reproduction)')
    rowsA = lll_rows(colsA, digits); genA = report(rowsA, 'A')
    colsB = [target]+kv+[r2*v for v in kv]+[Lchi, r2*Lchi]
    log('TEST B: target in extras-kernel span + L')
    rowsB = lll_rows(colsB, digits); genB = report(rowsB, 'B')
    colsC = [Lchi]+kv+[r2*v for v in kv]
    log('TEST C: L in extras-kernel span?')
    rowsC = lll_rows(colsC, digits); genC = report(rowsC, 'C')
    json.dump(dict(digits=digits,
                   A=[dict(vec=r['vec'],resid=mp.nstr(r['residual'],6)) for r in genA[:3]],
                   B=[dict(vec=r['vec'],resid=mp.nstr(r['residual'],6)) for r in genB[:3]],
                   C=[dict(vec=r['vec'],resid=mp.nstr(r['residual'],6)) for r in genC[:3]]),
              open(HERE/f'lchi_diag_{digits}.json','w'), indent=1)

if __name__=='__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 800)
