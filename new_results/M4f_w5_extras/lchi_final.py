"""M4e: final candidate = S0 particular + S0 sat kernel + L + sum a_z u_z, expanded to atoms.

Expansion: atom coeffs = particular_S0 + sum_j (a_j+s2 b_j) Ksat[j]
           + sum_z coeff_z * (Li5(z) - particular_z(S0))  [as atom contributions]
           + L coefficient (separate atom 'Lchi').
Saves full atom-level Q(sqrt2) coefficients to lchi_final_{digits}.json and prints height.
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
    sol = NS.solve(S, 'fin')
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
    dropk = {17, 27, 30}
    kindep = [j for j in range(35) if j not in dropk]
    zlist = [NS.EXTRA[i] for i in [0,2,3,4,5]]
    uzs=[]; parts=[]
    for z in zlist:
        t = W.li_symbol(5, z); xpz, xqz = solve_for_target(basis, cols, t)
        pv = mp.fsum((mp.mpf(x.numerator)/x.denominator + r2*mp.mpf(q.numerator)/q.denominator)*v
                     for x,q,v in zip(xpz,xqz,vals))
        zval = mp.mpf(evaluate(('li',5,z), digits+60).mid().str(digits+110, radius=False))
        uzs.append(zval - pv); parts.append((xpz,xqz))
    colv = [target] + [kv[j] for j in kindep] + [r2*kv[j] for j in kindep] + \
           [Lchi, r2*Lchi] + uzs + [r2*u for u in uzs]
    tol = mp.mpf(10)**(-digits-20)
    rows = lll_rows(colv, digits)
    gen = [r for r in rows if r['vec'][0]!=0 and abs(r['residual'])<tol]
    gen.sort(key=lambda r: r['maxcoef'])
    if not gen: log('NO relation'); return
    r = gen[0]; vec=r['vec']; c0=vec[0]
    log(f'c0={c0} maxcoef={r["maxcoef"]} resid={mp.nstr(r["residual"],6)}')
    nk = len(kindep); nu = len(zlist)
    # atom coefficients over (65 S0 atoms) + (5 Li5 extras) + Lchi
    ca = dict((nm, F(0)) for nm in names); cb = dict((nm, F(0)) for nm in names)
    for i,nm in enumerate(names): ca[nm]+=xp0[i]; cb[nm]+=xq0[i]
    for kk,j in enumerate(kindep):
        aj=F(-vec[1+kk],c0); bj=F(-vec[1+nk+kk],c0)
        for i,nm in enumerate(names):
            if K[j][i]: ca[nm]+=aj*K[j][i]; cb[nm]+=bj*K[j][i]
    aL=F(-vec[1+2*nk],c0); bL=F(-vec[1+2*nk+1],c0)
    enames=[]
    for k,z in enumerate(zlist):
        cz_a=F(-vec[1+2*nk+2+k],c0); cz_b=F(-vec[1+2*nk+2+nu+k],c0)
        nm=f'Li5({z})'; enames.append(nm)
        ca[nm]=ca.get(nm,F(0))+cz_a; cb[nm]=cb.get(nm,F(0))+cz_b
        xpz,xqz=parts[k]
        for i,n0 in enumerate(names):
            ca[n0]-=cz_a*xpz[i]; cb[n0]-=cz_a*xqz[i]
            ca[n0]-=cz_b*xpz[i]*0  # sqrt2 part handled via xqz
        # careful: cz = cz_a + s2 cz_b multiplies (Li5(z) - (xpz+s2 xqz).atoms)
        # contribution to atom n0: -cz*(xpz+s2 xqz) = -(cz_a*xpz - 2 cz_b*xqz) - s2(cz_a*xqz + cz_b*xpz)
        # redo properly below
    # redo extras contribution cleanly
    ca = dict((nm, F(0)) for nm in names); cb = dict((nm, F(0)) for nm in names)
    for i,nm in enumerate(names): ca[nm]+=xp0[i]; cb[nm]+=xq0[i]
    for kk,j in enumerate(kindep):
        aj=F(-vec[1+kk],c0); bj=F(-vec[1+nk+kk],c0)
        for i,nm in enumerate(names):
            if K[j][i]: ca[nm]+=aj*K[j][i]; cb[nm]+=bj*K[j][i]
    for k,z in enumerate(zlist):
        cz_a=F(-vec[1+2*nk+2+k],c0); cz_b=F(-vec[1+2*nk+2+nu+k],c0)
        nm=f'Li5({z})'; ca[nm]=ca.get(nm,F(0))+cz_a; cb[nm]=cb.get(nm,F(0))+cz_b
        xpz,xqz=parts[k]
        for i,n0 in enumerate(names):
            ca[n0]+= -(cz_a*xpz[i] - 2*cz_b*xqz[i])
            cb[n0]+= -(cz_a*xqz[i] + cz_b*xpz[i])
    alln = names + enames
    dens=[ca[n].denominator for n in alln]+[cb[n].denominator for n in alln]+[aL.denominator,bL.denominator]
    D=lcm(*dens); mx=max([abs(ca[n].numerator) for n in alln]+[abs(cb[n].numerator) for n in alln])
    log(f'atoms+extras={len(alln)} lcm-den={D} max|num|={mx}')
    log(f'L coefficient: {aL} + {bL} sqrt2')
    li5=[(n,ca[n],cb[n]) for n in alln if n.startswith('Li5') and (ca[n] or cb[n])]
    log(f'nonzero Li5 args: {len(li5)}')
    json.dump(dict(digits=digits,c0=c0,aL=str(aL),bL=str(bL),names=alln,
                   ca=[str(ca[n]) for n in alln], cb=[str(cb[n]) for n in alln],
                   resid=mp.nstr(r['residual'],8)),
              open(HERE/f'lchi_final_{digits}.json','w'), indent=1)
    log(f'saved lchi_final_{digits}.json')

if __name__=='__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 1200)
