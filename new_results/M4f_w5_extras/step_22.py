WORK='<project>/research/08_theory/M4f_w5_extras'
src = '''"""M4f shared: cached high-precision column values for the search.

Columns (over Q, Q(sqrt2) doubled as v and sqrt2*v):
  kv_j for j in kindep (saturated S0 symbol-kernel values), L(chi8,5),
  u_z = Li5(z) - p_z(S0 atoms) for each pool extra z.
All mpf at `digits+120` working precision. Cache: vals_cache_{digits}.pkl.
Precision is set AFTER all imports (trap).
"""
import json, math, pickle
from fractions import Fraction as F
from pathlib import Path
import sys
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
from flint import arb, ctx
import mpmath as mp
import w5_symbol as W
from w5_search import evaluate
import nice_search as NS

ZERO, ONE, Q, R, A, B = W.ZERO, W.ONE, W.Q, W.R, W.A, W.B

def set_prec(digits):
    ctx.prec = math.ceil((digits+140)*math.log2(10))
    mp.mp.dps = digits+120

def qnum(c):
    return mp.mpf(c.p.numerator)/c.p.denominator + mp.sqrt(2)*mp.mpf(c.q.numerator)/c.q.denominator

def extras_from_pool():
    """Pool elements that are NOT S0 Li5 arguments (the genuine extras), in pool order."""
    d = json.loads((HERE/'pool.json').read_text())
    s0 = set(NS.S0)
    out = []
    for e in d['pool']:
        x = Q(F(e['p']), F(e['q']))
        if x in s0: continue
        out.append(dict(x=x, xp=[F(c) for c in e['xp']], xq=[F(c) for c in e['xq']], h=e['h']))
    return out

def build(digits=1200):
    """Returns dict with names, vals (65 S0 atom values), K, kindep, kv, Lchi, target,
    extras (list of dicts with x, uz, xp, xq). Uses cache for atom values."""
    set_prec(digits)
    sol = NS.solve(list(NS.S0), 'build')
    names = sol['names']; basis = sol['basis']; xp0, xq0 = sol['xp'], sol['xq']
    ks = json.loads((HERE/'nice_ksat_s0.json').read_text()); K = ks['K']
    recipes = [(b[0], b[2]) for b in basis]
    cache = HERE/f'vals_cache_{digits}.pkl'
    if cache.exists():
        vals = pickle.loads(cache.read_bytes())
    else:
        vals = [mp.mpf(v.mid().str(digits+110, radius=False))
                for v in [evaluate(r, digits+60) for r in recipes]]
        cache.write_bytes(pickle.dumps(vals))
    r2 = mp.sqrt(2)
    Lchi = mp.mpf(json.loads((HERE/'Lchi8_5.json').read_text())['value'])
    cert = json.loads((HERE/'cobracket.json').read_text())
    fixedv = mp.mpf(0)
    for term in cert['correction']:
        x,y,c = (Q(*map(F,term[k])) for k in ['x','y','coefficient'])
        v = evaluate(('mpl',3,2,x,y), digits+60)
        fixedv += qnum(c)*mp.mpf(v.mid().str(digits+110, radius=False))
    saved = json.loads((HERE/'W5_5400.json').read_text())
    pval = mp.fsum(qnum(W.Q(x,q_))*v for x,q_,v in zip(xp0,xq0,vals))
    target = mp.mpf(saved['value']) - fixedv - pval
    kv = [mp.fsum(F(c)*v for c,v in zip(krow,vals)) for krow in K]
    dropk = {17,27,30}
    for j in dropk:
        assert abs(kv[j]) < mp.mpf(10)**(-digits-60), f'kv{j} not zero'
    kindep = [j for j in range(35) if j not in dropk]
    extras = []
    for e in extras_from_pool():
        z = e['x']
        pv = mp.fsum(qnum(W.Q(a,b_))*v for a,b_,v in zip(e['xp'],e['xq'],vals))
        zval = mp.mpf(evaluate(('li',5,z), digits+60).mid().str(digits+110, radius=False))
        extras.append(dict(x=z, uz=zval-pv, xp=e['xp'], xq=e['xq'], h=e['h']))
    return dict(names=names, basis=basis, vals=vals, K=K, kindep=kindep, kv=kv,
                Lchi=Lchi, target=target, extras=extras, xp0=xp0, xq0=xq0)

def lll_rows(vals, digits, label=''):
    import signal, time
    n=len(vals); C=mp.mpf(10)**digits
    rows=[[int(i==j) for j in range(n)]+[int(mp.nint(C*v))] for i,v in enumerate(vals)]
    signal.alarm(2700); t0=time.monotonic()
    red=flint.fmpz_mat(rows).lll(gram='exact'); signal.alarm(0)
    out=[]
    for i in range(n):
        vec=[int(red[i,j]) for j in range(n)]
        out.append(dict(vec=vec, maxcoef=max(map(abs,vec)),
                        residual=mp.fsum(c*v for c,v in zip(vec,vals))))
    print(f'  LLL {n}x{n}@{digits} {label}: {time.monotonic()-t0:.1f}s', flush=True)
    return out
'''
open(f'{WORK}/m4f_common.py','w').write(src)
print('wrote m4f_common.py')
