"""Values for W6-N6 stage: general Li_{s,t}(x,y) with rigorous tails, leaf cache, product values.
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
        v = json.loads((HERE/'Lchi8_5_2500.json').read_text())['value']
        assert len(v) >= digits+20, 'Lchi precision insufficient'
        return arb(v[:digits+20])
    raise ValueError(recipe)

def leaf_recipes():
    recs = [[('log', Q(2)), ('log', A)],
            [('pi2',), ('li', 2, B)],
            [('zeta', 3), ('li', 3, B), ('li', 3, -B), ('li', 3, R/2)],
            [('li', 4, z) for z in [B, -B, -B*B, R/2, Q(F(1,2)), ONE-R/2, 2-R]],
            [('li', 5, z) for z in [B, -B, -B*B, R/2, Q(F(1,2)), ONE-R/2, 2-R,
                                    Q(F(1,2))-R/2, Q(F(1,2))-R/4, Q(F(1,2))+R/4]]
            + [('zeta', 5), ('lchi',)]
            + [('mpl', 3, 2, x, y) for x, y in [(B, ONE), (R/2, ONE), (-B, ONE), (B, -ONE)]],
            [('li', 6, z) for z in BASE_LI6]]
    out = [r for g in recs for r in g]
    for m, n, x, y, c in n6_terms():
        out.append(('mpl', m, n, x, y))
    return out

def main(digits):
    ctx.prec = math.ceil((digits+140)*math.log2(10))
    mp.mp.dps = digits+120
    cache = HERE/f'w6_leaves_{digits}.pkl'
    recs = leaf_recipes()
    if cache.exists():
        saved = pickle.loads(cache.read_bytes())
        if isinstance(saved, dict):
            saved = list(saved.values())
        assert len(saved) == len(recs)
        leaves = {rec: mp.mpf(v) for rec, v in zip(recs, saved)}
        log(f'loaded leaf cache ({len(leaves)} leaves)')
    else:
        t0 = time.monotonic()
        leaves = {}
        for rec in recs:
            v = leaf_value(rec, digits+40)
            leaves[rec] = mp.mpf(v.mid().str(digits+110, radius=False))
            log(f'  leaf {rec} done ({time.monotonic()-t0:.1f}s)')
        cache.write_bytes(pickle.dumps([mp.nstr(leaves[r], digits+110) for r in recs]))
        log(f'leaf values computed ({time.monotonic()-t0:.1f}s)')
    # atom values for the full basis
    basis = atoms_w6(BASE_LI6)
    acache = HERE/f'w6_atomvals_{digits}.pkl'
    if acache.exists():
        vals = [mp.mpf(v) for v in pickle.loads(acache.read_bytes())]
        log(f'loaded atom value cache ({len(vals)})')
    else:
        t0 = time.monotonic()
        vals = []
        for name, _, recipe in basis:
            if recipe[0] == 'product':
                v = mp.mpf(1)
                for leaf in recipe[1]:
                    v *= leaves[leaf]
                vals.append(v)
            elif recipe[0] in ('li', 'mpl'):
                vals.append(leaves[recipe])
            else:
                raise ValueError(recipe)
        acache.write_bytes(pickle.dumps([mp.nstr(v, digits+110) for v in vals]))
        log(f'atom values computed ({time.monotonic()-t0:.1f}s)')
    return leaves, vals

if __name__ == '__main__':
    digits = int(sys.argv[1]) if len(sys.argv) > 1 else 1650
    leaves, vals = main(digits)
    ctx.prec = math.ceil((digits+140)*math.log2(10)); mp.mp.dps = digits+120
    def lv(rec):
        return mp.mpf(leaf_value(rec, digits+40).mid().str(digits+110, radius=False))
    d1 = leaves[('li', 6, R/2)] + lv(('li', 6, -R/2)) - leaves[('li', 6, Q(F(1,2)))]/32
    print('Li6 duplication control (1/sqrt2 pair vs 1/2):', mp.nstr(abs(d1), 5))
    d2 = leaves[('li', 6, B)] + leaves[('li', 6, -B)] - lv(('li', 6, B*B))/32
    print('Li6 duplication control (b pair vs b^2):', mp.nstr(abs(d2), 5))
    import w5_search as WS5
    ctx.prec = math.ceil((digits+140)*math.log2(10))
    ref = WS5.mpl(3, 2, B, ONE, digits)
    dd = leaves[('mpl', 3, 2, B, ONE)] - mp.mpf(ref.mid().str(digits+110, radius=False))
    print('Li32(b,1) mpl_general vs M4c mpl:', mp.nstr(abs(dd), 5))
    n6v = json.loads((HERE.parent/'M5a_next_coefficient/n6_values.json').read_text())
    mp.mp.dps = 120
    for key, rec in [('Li33_b1', ('mpl', 3, 3, B, ONE)), ('Li33_s1', ('mpl', 3, 3, R/2, ONE)),
                     ('Li33_b2', ('mpl', 3, 3, B, Q(2))), ('Li42_bs', ('mpl', 4, 2, B, R/2))]:
        print(key, 'vs M5a n6_values:', mp.nstr(abs(leaves[rec] - mp.mpf(n6v[key])), 5))
    mp.mp.dps = 80
    zz = Q(F(1,2))-R/4
    zf = mp.mpf(1)/2 - mp.sqrt(2)/4
    direct = mp.polylog(6, zf)
    print('Li6((2-sqrt2)/4) acb vs mpmath polylog:', mp.nstr(abs(leaves[('li', 6, zz)]-direct), 5))
