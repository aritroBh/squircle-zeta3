WORK='<project>/research/08_theory/M4f_w5_extras'
src = r'''"""M4f: inspect + numerically verify the best optimized form (full set = same point as (1..6)).
Check: value of the assembled atom vector vs W5_5400 at ~1250 digits (cached values).
Also print Li5/Lchi/zeta5 coefficients and term census.
"""
import json, pickle, sys
from fractions import Fraction as F
from math import lcm
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mpmath as mp
import m4f_common as MC
import nice_search as NS

def main():
    out = pickle.loads((HERE/'optimized2.pkl').read_bytes())
    Sfull = (0,1,2,3,4,5,6)
    sc, bv, tt = out[Sfull]
    print('score:', sc)
    sol = NS.solve(list(NS.S0), 'insp')
    names = sol['names']
    pool = json.loads((HERE/'pool.json').read_text())['pool']
    s0 = set(NS.S0)
    extras = [e for e in pool if NS.W.Q(F(e['p']),F(e['q'])) not in s0]
    # apply display reduction for the printed coefficients
    ib2 = names.index('Li5((3-2*sqrt2))'); ibp = names.index('Li5((-1+1*sqrt2))')
    ibm = names.index('Li5((1-1*sqrt2))'); im1 = names.index('Li5(-1)'); iz5 = names.index('z5')
    v = list(bv)
    for off in (0,65):
        a = v[ib2+off]
        if a: v[ibp+off] += 16*a; v[ibm+off] += 16*a; v[ib2+off] = F(0)
        a = v[im1+off]
        if a: v[iz5+off] += F(-15,16)*a; v[im1+off] = F(0)
    print('--- Li5 / zeta5 / Lchi coefficients (reduced) ---')
    for i,nm in enumerate(names):
        a,b = v[i], v[65+i]
        if (a or b) and (nm.startswith('Li5') or nm=='z5'):
            print(f'  {nm:26s} {a} + ({b}) s2')
    print(f'  {"Lchi8_5":26s} {v[130]} + ({v[131]}) s2')
    for k,e in enumerate(extras):
        a,b = v[132+k], v[139+k]
        if a or b:
            print(f'  Li5({e["p"]}+{e["q"]}s2){"":12s} {a} + ({b}) s2')
    # numerical verification at 1250 digits
    d = MC.build(1200)
    vals = d['vals']; r2 = mp.sqrt(2)
    total = mp.mpf(0)
    for i in range(65):
        a,b = bv[i], bv[65+i]
        if a or b:
            total += (mp.mpf(a.numerator)/a.denominator + r2*mp.mpf(b.numerator)/b.denominator)*vals[i]
    Lchi = d['Lchi']
    total += (mp.mpf(v[130].numerator)/v[130].denominator + r2*mp.mpf(v[131].numerator)/v[131].denominator)*Lchi
    # NOTE: use UNREDUCED bv for value; reduced v equals same value by exact identities, check both
    total2 = mp.mpf(0)
    for i in range(65):
        a,b = v[i], v[65+i]
        if a or b:
            total2 += (mp.mpf(a.numerator)/a.denominator + r2*mp.mpf(b.numerator)/b.denominator)*vals[i]
    for k,e in enumerate(extras):
        a,b = bv[132+k], bv[139+k]
        if a or b:
            zval = mp.mpf(MC.evaluate(('li',5,e['x']), 1260).mid().str(1300, radius=False))
            total += (mp.mpf(a.numerator)/a.denominator + r2*mp.mpf(b.numerator)/b.denominator)*zval
            total2 += (mp.mpf(a.numerator)/a.denominator + r2*mp.mpf(b.numerator)/b.denominator)*zval
    total2 += (mp.mpf(v[130].numerator)/v[130].denominator + r2*mp.mpf(v[131].numerator)/v[131].denominator)*Lchi
    cert = json.loads((HERE/'cobracket.json').read_text())
    n5v = mp.mpf(0)
    for term in cert['correction']:
        x,y,c = (NS.W.Q(*map(F,term[k])) for k in ['x','y','coefficient'])
        vv = MC.evaluate(('mpl',3,2,x,y), 1260)
        n5v += MC.qnum(c)*mp.mpf(vv.mid().str(1300, radius=False))
    saved = json.loads((HERE/'W5_5400.json').read_text())
    w5 = mp.mpf(saved['value'])
    r_unred = total + n5v - w5
    r_red = total2 + n5v - w5
    print(f'value check: unreduced - W5 = {mp.nstr(r_unred, 6)}')
    print(f'value check: reduced   - W5 = {mp.nstr(r_red, 6)}')
    # census
    nz = [(nm, bv[i], bv[65+i]) for i,nm in enumerate(names) if bv[i] or bv[65+i]]
    nprod = sum(1 for nm,_,_ in nz if '*' in nm)
    nli4 = sum(1 for nm,_,_ in nz if 'Li4' in nm)
    nw32 = sum(1 for nm,_,_ in nz if ('Li3' in nm or 'Li2' in nm or 'z3' in nm) and '*' in nm)
    print(f'atoms nonzero: {len(nz)} (products {nprod}, Li4* {nli4}, w3w2 {nw32})')
main()
'''
open(f'{WORK}/inspect_winner.py','w').write(src)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','inspect_winner.py'],
                   cwd=WORK, capture_output=True, text=True, timeout=900)
print(r.stdout)
print('STDERR', r.stderr[-1200:])
