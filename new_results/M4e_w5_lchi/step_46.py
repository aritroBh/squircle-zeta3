script = r'''"""M4e: test B = target in span(extras saturated kernel + L), expanded to atoms.

Uses nice_ksat.json (saturated kernel, S0+EXTRA[:6], 41 vectors, entries<=33).
Columns: [target_e] + kve_j, s2*kve_j (41) + [L, s2*L]. 85 columns.
Reports the W5 relation and expands to atom coefficients (particular_e + kernel combo + L).
Saves the expanded Q(sqrt2) coefficients per atom to lchi_testB_{digits}.json.
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

def main(digits=1200):
    ctx.prec = math.ceil((digits+140)*math.log2(10)); mp.mp.dps = digits+120
    S = list(NS.S0) + NS.EXTRA[:6]
    sol = NS.solve(S, 'testB')
    assert sol['ok']
    names = sol['names']; basis = sol['basis']; xp, xq = sol['xp'], sol['xq']
    ks = json.loads((HERE/'nice_ksat.json').read_text())
    assert ks['names'] == names
    K = ks['K']    # 41 saturated kernel vectors over 71 atoms
    log(f'atoms={len(names)} kernel={len(K)} max|Kentry|={max(max(map(abs,r)) for r in K)}')
    recipes = [(b[0], b[2]) for b in basis]
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
    cols = [target]+kv+[r2*v for v in kv]+[Lchi, r2*Lchi]
    tol = mp.mpf(10)**(-digits-20)
    rows = lll_rows(cols, digits)
    gen = [r for r in rows if r['vec'][0]!=0 and abs(r['residual'])<tol]
    gen.sort(key=lambda r: r['maxcoef'])
    if not gen:
        log(f'NO relation at {digits}'); return
    r = gen[0]
    vec = r['vec']; c0 = vec[0]
    log(f'W5 coeff c0={c0} maxcoef={r["maxcoef"]} resid={mp.nstr(r["residual"],6)}')
    log(f'  L coefs (last two of vec): {vec[-2]}, {vec[-1]}')
    # expand: target = -(1/c0)[ sum_{j} (vec[1+j] + s2 vec[1+41+j]) kv_j + (vec[83]+s2 vec[84]) L ]
    # W5 - N5 = particular + that.  Atom coeffs = particular + sum (a_j+s2 b_j) K[j]
    # coefficient (rational,sqrt2) per atom:
    ca = [F(0)]*len(names); cb = [F(0)]*len(names)
    for i,(x,q) in enumerate(zip(xp,xq)):
        ca[i]+=x; cb[i]+=q
    for j in range(41):
        aj = F(-vec[1+j], c0); bj = F(-vec[1+41+j], c0)
        for i in range(len(names)):
            if K[j][i]:
                ca[i]+=aj*K[j][i]; cb[i]+=bj*K[j][i]
    aL = F(-vec[83], c0); bL = F(-vec[84], c0)
    log(f'L(chi8,5) coefficient: {aL} + {bL} sqrt2')
    # report per-atom
    from math import lcm as L_
    dens = [c.denominator for c in ca+cb] + [aL.denominator, bL.denominator]
    log(f'lcm denom = {L_(*dens)}')
    mx = max(max(abs(c.numerator) for c in ca), max(abs(c.numerator) for c in cb))
    log(f'max |numerator| = {mx}')
    # which extras appear?
    for i, nm in enumerate(names):
        if nm.startswith('Li5') and (ca[i] or cb[i]):
            pass
    extra_names = [f'Li5({z})' for z in NS.EXTRA[:6]]
    used_extras = [nm for nm in extra_names if (ca[names.index(nm)] or cb[names.index(nm)])]
    log(f'extras with nonzero coeff: {used_extras}')
    data = dict(digits=digits, c0=c0, names=names,
                ca=[str(c) for c in ca], cb=[str(c) for c in cb],
                aL=str(aL), bL=str(bL), residual=mp.nstr(r['residual'],8))
    (HERE/f'lchi_testB_{digits}.json').write_text(json.dumps(data)+'\n')
    log(f'saved lchi_testB_{digits}.json')

if __name__=='__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 1200)
'''
M4E = '<project>/research/08_theory/M4e_w5_lchi'
open(M4E+'/lchi_testB.py','w').write(script)
import subprocess, time
PY = '<project>/.venv/bin/python'
t0=time.monotonic()
r = subprocess.run(['nice', PY, 'lchi_testB.py', '1200'], cwd=M4E, capture_output=True, text=True, timeout=2700)
print(r.stdout); print("STDERR:", r.stderr[-1500:]); print(f"elapsed {time.monotonic()-t0:.0f}s")
open(M4E+'/progress.log','a').write("step_44: lchi_testB.py 1200 (test B + atom expansion)\n")
