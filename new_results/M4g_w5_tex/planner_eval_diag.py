"""Planner check of M4e's R5_nice, with independent mpmath code (no imports from the M4c/M4d/M4e scripts).

R5_nice = N5 + sum (ca + cb sqrt2) atom + (aL + bL sqrt2) L(chi8,5), from lchi_final_{1200,1500}.json.
- Li_k(x), |x| < 1: direct power series. Li5(-1) = -(15/16) zeta(5).
- Li_{3,2}(x,y) = sum_{n>m>0} x^n y^m/(n^3 m^2): direct series with a running inner sum.
- N5: hard-coded from the M4d/M4e briefs (not read from cobracket.json).
- L(chi8,5) = 8^-5 [zeta(5,1/8) + zeta(5,7/8) - zeta(5,3/8) - zeta(5,5/8)] (mpmath Hurwitz), compared with
  M4e's Lchi8_5.json.
Also builds the reduced representative: Li5(b^2) -> 16 Li5(b) + 16 Li5(-b), Li5(-1) -> -(15/16) zeta(5).
Compares with W5_5400.json (M4c, certified arb) and W5_600.txt (M4a, a different quadrature code).
Usage: python planner_eval_R5nice.py [DPS]
"""
import json, re, sys
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp

HERE = Path('../M4e_w5_lchi').resolve()
TH = Path('..').resolve()
DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 1600
mp.mp.dps = DPS + 40
r2 = mp.sqrt(2)
lam, L = mp.log(2), mp.log(1 + r2)
EPS = mp.mpf(10) ** (-(DPS + 30))


def field(s):  # '(-1+1*sqrt2)', '1/2*sqrt2', '-1', '(1/2-3/8*sqrt2)'
    s = s.strip().strip('()').replace('sqrt2', 's')
    v = mp.mpf(0)
    for t in re.findall(r'[+-]?[^+-]+', s):
        if 's' in t:
            c = t.replace('*s', '').replace('s', '')
            c = c if c not in ('', '+', '-') else c + '1'
            v += mp.mpf(F(c).numerator) / F(c).denominator * r2
        else:
            v += mp.mpf(F(t).numerator) / F(t).denominator
    return v


def li(k, x):
    if x == -1:
        return -(1 - mp.mpf(2) ** (1 - k)) * mp.zeta(k)
    assert abs(x) < 1
    tot, xn, n = mp.mpf(0), mp.mpf(1), 1
    while True:
        xn *= x
        tot += xn / mp.mpf(n) ** k
        if abs(xn) < EPS: return tot
        n += 1


def li32(x, y):
    tot, H, xn, ym, n = mp.mpf(0), mp.mpf(0), mp.mpf(1), mp.mpf(1), 1
    while True:
        xn *= x
        tot += xn * H / mp.mpf(n) ** 3
        if n > 50 and abs(xn) < EPS: return tot
        ym *= y
        H += ym / mp.mpf(n) ** 2
        n += 1


cache = {}
def atom(a):
    if a not in cache:
        if a == 'lam': v = lam
        elif a == 'L': v = L
        elif a == 'pi2': v = mp.pi ** 2
        elif a == 'z3': v = mp.zeta(3)
        elif a == 'z5': v = mp.zeta(5)
        elif a == 'Li2(b)': v = li(2, r2 - 1)
        else:
            m = re.fullmatch(r'Li([2-5])\((.*)\)', a)
            assert m, a
            v = li(int(m.group(1)), field(m.group(2)))
        cache[a] = v
    return cache[a]


def split(name):
    out, depth, cur = [], 0, ''
    for ch in name:
        depth += (ch == '(') - (ch == ')')
        if ch == '*' and depth == 0: out.append(cur); cur = ''
        else: cur += ch
    return out + [cur]


def num(p, q): return mp.mpf(p.numerator) / p.denominator + r2 * mp.mpf(q.numerator) / q.denominator

b = r2 - 1
N5 = (num(F(779), F(-695, 56)) * li32(b, 1) - 672 * li32(1 / r2, 1)
      - num(F(374), F(-76)) * li32(-b, 1) - num(F(405), F(-4951, 56)) * li32(b, -1))

h = [mp.zeta(5, mp.mpf(k) / 8) for k in (1, 3, 5, 7)]
Lchi = (h[0] + h[3] - h[1] - h[2]) / mp.mpf(8) ** 5
Lsaved = mp.mpf(json.loads((HERE / 'Lchi8_5.json').read_text())['value'])
print(f'dps {DPS}')
print('L(chi8,5) =', mp.nstr(Lchi, 40))
print('L(chi8,5) - M4e Lchi8_5.json =', mp.nstr(Lchi - Lsaved, 6))

W5c = mp.mpf(json.loads((HERE / 'W5_5400.json').read_text())['value'])
W5a = mp.mpf(re.sub(r'\s', '', (TH / 'M4a_decomposition' / 'W5_600.txt').read_text()).split('=')[-1])


def reduce_rep(names, ca, cb):
    """Li5(b^2) = 16 Li5(b) + 16 Li5(-b); Li5(-1) = -(15/16) zeta(5)."""
    c = {nm: [a, q] for nm, a, q in zip(names, ca, cb)}
    B2, Bp, Bm = 'Li5((3-2*sqrt2))', 'Li5((-1+1*sqrt2))', 'Li5((1-1*sqrt2))'
    a, q = c.pop(B2)
    c[Bp][0] += 16 * a; c[Bp][1] += 16 * q; c[Bm][0] += 16 * a; c[Bm][1] += 16 * q
    a, q = c.pop('Li5(-1)')
    c['z5'][0] += F(-15, 16) * a; c['z5'][1] += F(-15, 16) * q
    return {nm: v for nm, v in c.items() if v[0] or v[1]}


reps = {}
for D in (1200, 1500):
    d = json.loads((HERE / f'lchi_final_{D}.json').read_text())
    names, ca, cb = d['names'], [F(x) for x in d['ca']], [F(x) for x in d['cb']]
    aL, bL = F(d['aL']), F(d['bL'])
    for tag, coeffs in (('as given', {nm: (a, q) for nm, a, q in zip(names, ca, cb) if a or q}),
                        ('reduced', reduce_rep(names, ca, cb))):
        tot = N5 + num(aL, bL) * Lchi
        for nm, (a, q) in coeffs.items():
            v = mp.mpf(1)
            for f in split(nm): v *= atom(f)
            tot += num(a, q) * v
        reps[(D, tag)] = tot
        nli5 = sum(1 for nm in coeffs if nm.startswith('Li5('))
        dens = [x.denominator for v in coeffs.values() for x in v] + [aL.denominator, bL.denominator]
        nums = [abs(x.numerator) for v in coeffs.values() for x in v] + [abs(aL.numerator), abs(bL.numerator)]
        print(f'[{D} {tag}] terms {len(coeffs) + 5} (Li5 {nli5}, z5 {"z5" in coeffs}), '
              f'max |num| {max(nums)}, lcm den {mp.mpf(__import__("math").lcm(*dens))}')
        print(f'   R5_nice - W5_M4c(5400 bits) = {mp.nstr(tot - W5c, 6)}')
        print(f'   R5_nice - W5_M4a(600 digits) = {mp.nstr(tot - W5a, 6)}')
print('rep1200 - rep1500 =', mp.nstr(reps[(1200, 'as given')] - reps[(1500, 'as given')], 6))

# --- added block: evaluate planner/R5_nice_reduced.json directly ---
d2 = json.loads((HERE / 'planner' / 'R5_nice_reduced.json').read_text())['terms']
tot2 = N5
for nm, (a, q) in d2.items():
    if nm == 'Lchi8_5':
        v = Lchi
    else:
        v = mp.mpf(1)
        for f in split(nm):
            v *= atom(f)
    tot2 += num(F(a), F(q)) * v
print('saved R5_nice_reduced.json - W5_M4c =', mp.nstr(tot2 - W5c, 6))
print('N5 =', mp.nstr(N5, 30))
for f in ('Li5((-1+1*sqrt2))', 'Li5((1/2+1/4*sqrt2))'):
    print(f, mp.nstr(atom(f), 30))

# dump all needed atom values to atoms_planner.json
need = set()
for nm in d2:
    if nm == 'Lchi8_5': continue
    for f in split(nm): need.add(f)
out = {}
for f in sorted(need):
    v = atom(f)
    out[f] = mp.nstr(v, 460)
out['Lchi8_5'] = mp.nstr(Lchi, 460)
out['N5'] = mp.nstr(N5, 460)
json.dump(out, open('atoms_planner.json','w'))
print('dumped', len(out))
