script = r'''"""M4e: identify the 2nd missing weight-5 constant.

For each extras Li5 argument z: solve symbol(Li5(z)) over the S0 atoms
(particular p_z), form the symbol-zero value u_z = Li5(z) - val(p_z),
and test u_z for membership in span(saturated S0 kernel values, L, s2*L)
via LLL. Also print the exact saturated-kernel zero-value relations.
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
    return out

def solve_for_target(basis, cols, t, label):
    words = NS.WORDS
    tp = [F(str(t.get(w, ZERO).p)) for w in words]
    tq = [F(str(t.get(w, ZERO).q)) for w in words]
    rows = [[str(c) for c in r] + [str(tp[i]), str(tq[i])] for i, r in enumerate(zip(*cols))]
    mat, rank = flint.fmpq_mat(rows).rref()
    n = len(basis)
    pivots = [next(j for j in range(mat.ncols()) if mat[i, j]) for i in range(rank)]
    ok = all(j < n for j in pivots)
    if not ok:
        return None
    xp = [F(0)]*n; xq = [F(0)]*n
    for i, p in enumerate(pivots):
        xp[p] = F(str(mat[i, n])); xq[p] = F(str(mat[i, n+1]))
    return xp, xq

def main(digits=800):
    ctx.prec = math.ceil((digits+140)*math.log2(10)); mp.mp.dps = digits+120
    S = list(NS.S0)
    sol = NS.solve(S, 'uz')
    names = sol['names']; basis = sol['basis']
    cols = NS.symbol_columns(basis)
    ks = json.loads((HERE/'nice_ksat_s0.json').read_text()); K = ks['K']
    recipes = [(b[0], b[2]) for b in basis]
    values = [evaluate(r, digits+60) for _, r in recipes]
    vals = [mp.mpf(v.mid().str(digits+110, radius=False)) for v in values]
    r2 = mp.sqrt(2)
    Lchi = mp.mpf(json.loads((HERE/'Lchi8_5.json').read_text())['value'])
    kv = [mp.fsum(F(c)*v for c,v in zip(krow,vals)) for krow in K]
    # exact zero-value kernel relations
    log('exact zero-value saturated-kernel vectors (resid < 1e-%d):' % (digits+80))
    for j, v in enumerate(kv):
        if abs(v) < mp.mpf(10)**(-digits-80):
            terms = [(names[i], c) for i, c in enumerate(K[j]) if c]
            log(f'  kv{j} = 0 exactly: {terms}')
    # independent V0 basis: drop kv30,s2kv30 (exact 0); keep rest, LLL handles remaining deps
    drop = {30}
    indep = [j for j in range(35) if j not in drop]
    basecols = [kv[j] for j in indep] + [r2*kv[j] for j in indep] + [Lchi, r2*Lchi]
    tol = mp.mpf(10)**(-digits-20)
    for z in NS.EXTRA[:6]:
        t = W.li_symbol(5, z)
        out = solve_for_target(basis, cols, t, str(z))
        if out is None:
            log(f'z={z}: symbol NOT in S0 span'); continue
        xp, xq = out
        pv = mp.fsum((mp.mpf(x.numerator)/x.denominator + r2*mp.mpf(q.numerator)/q.denominator)*v
                     for x, q, v in zip(xp, xq, vals))
        zval = mp.mpf(evaluate(('li', 5, z), digits+60).mid().str(digits+110, radius=False))
        uz = zval - pv
        cols_ = [uz] + basecols
        rows = lll_rows(cols_, digits)
        gen = [r for r in rows if r['vec'][0] != 0 and abs(r['residual']) < tol]
        gen.sort(key=lambda r: r['maxcoef'])
        status = (f'IN SPAN: coeff {gen[0]["vec"][0]}, maxcoef {gen[0]["maxcoef"]}' if gen
                  else 'NOT in span (new constant direction)')
        log(f'z={z}: u_z ~ {mp.nstr(uz, 10)} -> {status}')
        if gen:
            r = gen[0]
            log(f'    vec[:8]={r["vec"][:8]} Lcoefs={r["vec"][-2:]} resid={mp.nstr(r["residual"],4)}')

if __name__=='__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 800)
'''
M4E = '<project>/research/08_theory/M4e_w5_lchi'
open(M4E+'/lchi_second.py','w').write(script)
import subprocess, time
PY = '<project>/.venv/bin/python'
t0=time.monotonic()
r = subprocess.run(['nice', PY, 'lchi_second.py', '800'], cwd=M4E, capture_output=True, text=True, timeout=2700)
print(r.stdout); print("STDERR:", r.stderr[-1500:]); print(f"elapsed {time.monotonic()-t0:.0f}s")
open(M4E+'/progress.log','a').write("step_43: lchi_second.py 800 (identify 2nd constant)\n")
