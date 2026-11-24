import os
m5b='<project>/research/08_theory/M5b_w6_closed_form'
script = '''"""Values for W6-N6 stage: general Li_{s,t}(x,y) with rigorous tails, leaf cache, product values.
Usage: -B w6_values.py <digits>   (computes leaf + atom values, caches pkl)
"""
import json, math, pickle, sys, time
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
from flint import arb, acb, ctx
import mpmath as mp
from w4_symbol import Q, R, A, B
from m5b_common import ZERO, ONE, log
from w6_solve import atoms_w6, n6_terms, BASE_LI6

def numeric(z):
    return arb(z.p.numerator)/z.p.denominator + arb(2).sqrt()*arb(z.q.numerator)/z.q.denominator

def mpl_general(s, t, x, y, digits):
    """Li_{s,t}(x,y) = sum_{n>m>0} x^n y^m/(n^s m^t), for |x|<1, |x*y|<1. Rigorous ball tail."""
    xx, yy = numeric(x), numeric(y)
    rx, ry = abs(xx), abs(yy)
    assert 0 < rx < 1 and rx*ry < 1
    r = rx*ry
    N = int((digits+30)*math.log(10)/-math.log(float(r))) + 10
    total, H, xp, yp = arb(0), arb(0), arb(1), arb(1)
    for n in range(1, N+1):
        xp *= xx
        total += xp*H/(n**s)
        yp *= yy
        H += yp/(n**t)
    # tail: |H(n)| <= y^n/(y-1) if y>1 else n ; tail <= sum_{n>N} rx^n/n^s |H(n)|
    if ry > 1:
        tail = r**(N+1)/((ry-1)*(N+1)**s*(1-r))
    else:
        tail = rx**(N+1)/((N+1)**(s-1)*(1-rx))
    return total + arb(0, tail.abs_upper())

def leaf_value(recipe, digits):
    kind, *a = recipe
    if kind == 'log':
        return numeric(a[0]).log()
    if kind == 'pi2':
        return arb.pi()*arb.pi()
    if kind == 'zeta':
        return arb(a[0]).zeta()
    if kind == 'li':
        z = acb(numeric(a[1])).polylog(a[0])
        assert z.imag.contains(0)
        return z.real
    if kind == 'mpl':
        return mpl_general(a[0], a[1], a[2], a[3], digits)
    if kind == 'lchi':
        v = json.loads((HERE/'Lchi8_5.json').read_text())['value']
        assert len(v) >= digits, 'Lchi8_5.json precision insufficient'
        return arb(v[:digits+20])
    raise ValueError(recipe)

def get_leaves(digits):
    leaves = {}
    def add(rec):
        if rec not in leaves:
            leaves[rec] = leaf_value(rec, digits)
    groups = [
        [('log', Q(2)), ('log', A)],
        [('pi2',), ('li', 2, B)],
        [('zeta', 3), ('li', 3, B), ('li', 3, -B), ('li', 3, R/2)],
        [('li', 4, z) for z in [B, -B, -B*B, R/2, Q(F(1,2)), ONE-R/2, 2-R]],
        [('li', 5, z) for z in [B, -B, -B*B, R/2, Q(F(1,2)), ONE-R/2, 2-R,
                                Q(F(1,2))-R/2, Q(F(1,2))-R/4, Q(F(1,2))+R/4]]
        + [('zeta', 5), ('lchi',)]
        + [('mpl', 3, 2, x, y) for x, y in [(B, ONE), (R/2, ONE), (-B, ONE), (B, -ONE)]],
        [('li', 6, z) for z in BASE_LI6],
    ]
    for g in groups:
        for rec in g:
            add(rec)
    # N6 atoms
    for m, n, x, y, c in n6_terms():
        add(('mpl', m, n, x, y))
    return leaves

def main(digits):
    ctx.prec = math.ceil((digits+140)*math.log2(10))
    mp.mp.dps = digits+120
    cache = HERE/f'w6_leaves_{digits}.pkl'
    if cache.exists():
        leaves_str = pickle.loads(cache.read_bytes())
        leaves = {eval(k): v for k, v in leaves_str.items()}
        log(f'loaded leaf cache ({len(leaves)} leaves)')
    else:
        t0 = time.monotonic()
        leaves = {}
        for rec in list(get_leaves(digits)):
            leaves[rec] = leaf_value(rec, digits+40)
            log(f'  leaf {rec} done ({time.monotonic()-t0:.1f}s)')
        cache.write_bytes(pickle.dumps({repr(k): v for k, v in leaves.items()}))
        log(f'leaf values computed ({time.monotonic()-t0:.1f}s)')
    # atom values for the full basis
    basis = atoms_w6(BASE_LI6)
    acache = HERE/f'w6_atomvals_{digits}.pkl'
    if acache.exists():
        vals = pickle.loads(acache.read_bytes())
        log(f'loaded atom value cache ({len(vals)})')
    else:
        t0 = time.monotonic()
        vals = []
        for name, _, recipe in basis:
            if recipe[0] == 'product':
                v = mp.mpf(1)
                for leaf in recipe[1]:
                    v *= mp.mpf(leaves[leaf].mid().str(digits+110, radius=False))
                vals.append(v)
            elif recipe[0] == 'li':
                vals.append(mp.mpf(leaves[recipe].mid().str(digits+110, radius=False)))
            else:
                raise ValueError(recipe)
        acache.write_bytes(pickle.dumps(vals))
        log(f'atom values computed ({time.monotonic()-t0:.1f}s)')
    return leaves, vals

if __name__ == '__main__':
    digits = int(sys.argv[1]) if len(sys.argv) > 1 else 1650
    leaves, vals = main(digits)
    # controls: duplication Li6(z)+Li6(-z)=Li6(z^2)/32 for z=1/sqrt2 (z^2=1/2)
    z = R/2
    lhs = leaves[('li', 6, z)] + leaves[('li', 6, -z)] if ('li', 6, -z) in leaves else None
    d1 = leaves[('li', 6, R/2)] + leaves[('li', 6, -R/2)] - leaves[('li', 6, Q(F(1,2)))]/32
    print('Li6 duplication control (1/sqrt2 pair vs 1/2):', d1.str(8))
    d2 = leaves[('li', 6, B)] + leaves[('li', 6, -B)] - leaves[('li', 6, B*B)]/32
    print('Li6 duplication control (b pair vs b^2):', d2.str(8))
    d3 = leaves[('mpl', 1, 1, B, ONE)] if ('mpl',1,1,B,ONE) in leaves else None
    # Li32(b,1) vs M4c mpl
    import w5_search as WS5
    ref = WS5.mpl(3, 2, B, ONE, digits)
    dd = leaves[('mpl', 3, 2, B, ONE)] - ref
    print('Li32(b,1) vs M4c mpl:', dd.str(8))
    # N6 atoms vs M5a n6_values.json (90 digits)
    n6v = json.loads((HERE.parent/'M5a_next_coefficient/n6_values.json').read_text())
    mp.mp.dps = 120
    for key, rec in [('Li33_b1', ('mpl', 3, 3, B, ONE)), ('Li33_s1', ('mpl', 3, 3, R/2, ONE)),
                     ('Li33_b2', ('mpl', 3, 3, B, Q(2))), ('Li42_bs', ('mpl', 4, 2, B, R/2))]:
        got = mp.mpf(leaves[rec].mid().str(110, radius=False))
        print(key, 'vs M5a:', mp.nstr(abs(got - mp.mpf(n6v[key])), 5))
    # Li6 vs acb direct cross-check at z=1/2-1/4 sqrt2 via mpmath series
    mp.mp.dps = 80
    zz = Q(F(1,2))-R/4
    zf = mp.mpf(str(zz.p.numerator))/zz.p.denominator + mp.sqrt(2)*mp.mpf(str(zz.q.numerator))/zz.q.denominator
    direct = mp.polylog(6, zf)
    got = mp.mpf(leaves[('li', 6, zz)].mid().str(70, radius=False))
    print('Li6((2-sqrt2)/4) acb vs mpmath polylog:', mp.nstr(abs(got-direct), 5))
'''
open(os.path.join(m5b,'w6_values.py'),'w').write(script)
import py_compile; py_compile.compile(os.path.join(m5b,'w6_values.py'), doraise=True)
print('ok')
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'-B','w6_values.py','1650'],capture_output=True,text=True,cwd=m5b,timeout=1700)
print(r.stdout[-3000:]); print(r.stderr[-900:])
