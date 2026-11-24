"""Independent value check (Task 4): evaluate the final W5'' and W6 forms with direct-series
mpmath code only (no shared evaluator; no flint), at >=300 digits, against the quadrature values.
Li_s(z) = sum z^n/n^s (|z|<1); Li_{s,t}(x,y) = sum_{n>m>0} x^n y^m/(n^s m^t).
"""
import json, math, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mpmath as mp

DPS = 330
mp.mp.dps = DPS + 60
R2 = mp.sqrt(2)

def q2p(s):
    """parse QS2 repr like (-1+1*sqrt2) or 1/2*sqrt2 to mpf."""
    s = s.strip().strip('()').replace(' ', '')
    if 'sqrt2' not in s:
        return mp.mpf(s)
    s = s.replace('*sqrt2', 'S')
    out = mp.mpf(0)
    import re
    for t in re.findall(r'[+-]?[^+-]+', s):
        if t.endswith('S'):
            c = t[:-1]
            if c in ('', '+'): c = '1'
            elif c == '-': c = '-1'
            out += mp.mpf(c)*R2
        else:
            out += mp.mpf(t)
    return out

def ds_li(s, z):
    z = mp.mpf(z)
    az = abs(z)
    assert az < 1
    N = int((DPS+80)*mp.log(10)/-mp.log(az)) + 20
    total = mp.mpf(0)
    power = mp.mpf(1)
    for n in range(1, N+1):
        power *= z
        total += power/n**s
    return total

def ds_mpl(s, t, x, y):
    x, y = mp.mpf(x), mp.mpf(y)
    assert abs(x) < 1 and abs(x*y) < 1
    r = abs(x*y)
    N = int((DPS+80)*mp.log(10)/-mp.log(r)) + 20
    total, H, xp, yp = mp.mpf(0), mp.mpf(0), mp.mpf(1), mp.mpf(1)
    for n in range(1, N+1):
        xp *= x
        total += xp*H/n**s
        yp *= y
        H += yp/n**t
    return total

def lchi5():
    # project convention: chi of Q(sqrt2), +1 on n = +/-1 mod 8, -1 on n = +/-3 mod 8
    return (mp.zeta(5, mp.mpf(1)/8) - mp.zeta(5, mp.mpf(3)/8)
            - mp.zeta(5, mp.mpf(5)/8) + mp.zeta(5, mp.mpf(7)/8))/8**5

def leaf_value(leaf):
    kind = leaf[0]
    if kind == 'log':
        return mp.log(q2p(leaf[1]))
    if kind == 'pi2':
        return mp.pi**2
    if kind == 'zeta':
        return mp.zeta(leaf[1])
    if kind == 'li':
        return ds_li(leaf[1], q2p(leaf[2]))
    if kind == 'mpl':
        return ds_mpl(leaf[1], leaf[2], q2p(leaf[3]), q2p(leaf[4]))
    if kind == 'lchi':
        return lchi5()
    raise ValueError(leaf)

def name_to_recipes(names_basis):
    """map atom name -> recipe, using the same JSON-stored string forms the files use."""
    out = {}
    for nm, rec in names_basis:
        out[nm] = rec
    return out

def ser_rec(rec):
    """recipe with QS2 objects -> recipe with strings (for leaf_value)."""
    def sv(v):
        return str(v) if hasattr(v, 'p') else v
    return tuple([rec[0]] + [sv(v) for v in rec[1:]])

def evaluate_form(cf, name2rec):
    total = mp.mpf(0)
    for nm, a, b in zip(cf['names'], cf['ca'], cf['cb']):
        coeff = mp.mpf(a) + R2*mp.mpf(b)
        if not coeff:
            continue
        rec = name2rec[nm]
        if rec[0] == 'product':
            v = mp.mpf(1)
            for leaf in rec[1]:
                v *= leaf_value(ser_rec(leaf))
        else:
            v = leaf_value(ser_rec(rec))
        total += coeff*v
    return total

def main():
    # ---- W5'': R5'' + N5'' vs quadrature W5pp_4800
    from m5b_common import pool_zs, atoms_w5
    basis5 = atoms_w5(pool_zs())
    n2r5 = {b[0]: b[2] for b in basis5}
    cf5 = json.loads((HERE/'w5pp_closed_form.json').read_text())
    tot5 = evaluate_form(cf5, n2r5)
    n5 = json.loads((HERE/'cobracket5pp.json').read_text())
    for t in n5['correction']:
        coeff = mp.mpf(t['coefficient'][0]) + R2*mp.mpf(t['coefficient'][1])
        xv = mp.mpf(t['x'][0]) + R2*mp.mpf(t['x'][1])
        yv = mp.mpf(t['y'][0]) + R2*mp.mpf(t['y'][1])
        tot5 += coeff*ds_mpl(3, 2, xv, yv)
    w5pp = mp.mpf(json.loads((HERE/'w5pp_4800.json').read_text())['value_full'])
    r5 = tot5 - w5pp
    print(f"W5'': direct-series R5''+N5'' - quadrature = {mp.nstr(r5, 6)} (300+ digits check)")
    # ---- W6: R6 + N6 vs quadrature w6_5800
    from w6_solve import atoms_w6, n6_terms, BASE_LI6
    basis6 = atoms_w6(BASE_LI6)
    n2r6 = {b[0]: b[2] for b in basis6}
    cf6 = json.loads((HERE/'w6_closed_form.json').read_text())
    tot6 = evaluate_form(cf6, n2r6)
    for m, n, x, y, c in n6_terms():
        coeff = mp.mpf(str(c.p)) + R2*mp.mpf(str(c.q))
        tot6 += coeff*ds_mpl(m, n, q2p(str(x)), q2p(str(y)))
    w6 = mp.mpf(json.loads((HERE/'w6_5800.json').read_text())['value_full'])
    r6 = tot6 - w6
    print(f'W6: direct-series R6+N6 - quadrature = {mp.nstr(r6, 6)} (300+ digits check)')
    print('lchi5 direct:', mp.nstr(lchi5(), 30))

if __name__ == '__main__':
    main()
