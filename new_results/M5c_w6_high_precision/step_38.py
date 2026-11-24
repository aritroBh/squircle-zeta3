indep='''"""M5c independent direct-series check for W6 = R6 + N6 (no shared evaluator, mpmath only).
Li_s(z) = sum z^n/n^s; Li_{s,t}(x,y) = sum_{n>m>0} x^n y^m/(n^s m^t); lchi via Hurwitz zeta.
Usage: -B w6c_indep_check.py [tag]
"""
import json, math, sys, re
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mpmath as mp

DPS = 330
mp.mp.dps = DPS + 80
R2 = mp.sqrt(2)

def q2p(s):
    s = s.strip().strip('()').replace(' ', '')
    if 'sqrt2' not in s:
        return mp.mpf(s)
    s = s.replace('*sqrt2', 'S')
    out = mp.mpf(0)
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
    z = mp.mpf(z); az = abs(z)
    assert az < 1
    N = int((DPS+80)*mp.log(10)/-mp.log(az)) + 20
    total = mp.mpf(0); power = mp.mpf(1)
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
    return (mp.zeta(5, mp.mpf(1)/8) - mp.zeta(5, mp.mpf(3)/8)
            - mp.zeta(5, mp.mpf(5)/8) + mp.zeta(5, mp.mpf(7)/8))/8**5

def leaf_value(leaf):
    kind = leaf[0]
    if kind == 'log': return mp.log(q2p(leaf[1]))
    if kind == 'pi2': return mp.pi**2
    if kind == 'zeta': return mp.zeta(leaf[1])
    if kind == 'li': return ds_li(leaf[1], q2p(leaf[2]))
    if kind == 'mpl': return ds_mpl(leaf[1], leaf[2], q2p(leaf[3]), q2p(leaf[4]))
    if kind == 'lchi': return lchi5()
    raise ValueError(leaf)

def ser_rec(rec):
    def sv(v):
        return str(v) if hasattr(v, 'p') else v
    return tuple([rec[0]] + [sv(v) for v in rec[1:]])

def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else ''
    from w6_solve import atoms_w6, n6_terms, BASE_LI6, pool_zs
    li6 = list(BASE_LI6)
    if tag:
        pool = pool_zs()
        for i in [int(i) for i in tag[3:].split('_')]:
            if pool[i] not in li6:
                li6.append(pool[i])
    basis6 = atoms_w6(li6)
    n2r6 = {b[0]: b[2] for b in basis6}
    cf6 = json.loads((HERE/'w6_closed_form.json').read_text())
    total = mp.mpf(0)
    for nm, a, b in zip(cf6['names'], cf6['ca'], cf6['cb']):
        coeff = mp.mpf(a) + R2*mp.mpf(b)
        if not coeff:
            continue
        rec = n2r6[nm]
        if rec[0] == 'product':
            v = mp.mpf(1)
            for leaf in rec[1]:
                v *= leaf_value(ser_rec(leaf))
        else:
            v = leaf_value(ser_rec(rec))
        total += coeff*v
    for m, n, x, y, c in n6_terms():
        coeff = mp.mpf(str(c.p)) + R2*mp.mpf(str(c.q))
        total += coeff*ds_mpl(m, n, q2p(str(x)), q2p(str(y)))
    w6 = mp.mpf(json.loads((HERE/'w6_11800.json').read_text())['value_full'])
    r6 = total - w6
    print(f'W6: direct-series R6+N6 - quadrature(11800 bits) = {mp.nstr(r6, 6)} (at {DPS} digits, no shared evaluator)')

if __name__ == '__main__':
    main()
'''
open('w6c_indep_check.py','w').write(indep)
import py_compile; py_compile.compile('w6c_indep_check.py', doraise=True)
print('w6c_indep_check.py ok')
import subprocess, os
work='<project>/research/08_theory/M5c_w6_high_precision'
print(subprocess.run(['bash','-lc',f'tail -6 "{work}/w6_vals_3500.out"'],capture_output=True,text=True).stdout)
try: os.kill(int(open(work+'/w6_vals_3500.pid').read()),0); print('vals job alive')
except OSError: print('vals job done')
