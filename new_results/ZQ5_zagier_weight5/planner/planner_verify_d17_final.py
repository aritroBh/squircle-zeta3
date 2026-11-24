# Planner verification of ZQ5b's minimised D=17 form (ZQ5b_zagier_w5_larger_S/d17_final.json):
# (1) value: T_coef*L(chi17,5)/sqrt17 + z5_coef*zeta(5) + sum c*Li5(x) + sum c*products = 0, with own Li5 (planner_li5),
#     own Hurwitz L-value and own letter logs, at two precisions;
# (2) S-unit factorisations of x and 1-x over the letters, found here by PSLQ (not read from ZQ5b);
# (3) exact symbol of the right-hand side in (Z^4)^{(x)5} is ZERO (pi^2, zeta(5), T have zero symbol);
# (4) stats: Li5 terms, product terms, max numerator, lcm denominator.
import sys, json, itertools
from fractions import Fraction as F
from collections import defaultdict
from math import lcm
import mpmath as mp
sys.path.insert(0, '<project>/research/08_theory/ZQ5_zagier_weight5/planner')
from planner_li5 import li5

J = json.load(open('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S/d17_final.json'))
LET = ['lE17', 'ls17', 'lp', "lp'"]


def parse(s):  # 'a+b*sqrt17' with a, b rational strings
    s = s.replace('+-', '-')
    i = s.rfind('*sqrt17')
    head = s[:i]
    k = max(head.rfind('+'), head.rfind('-'))
    return F(head[:k]), F(head[k:])


args = [(parse(x), F(c)) for x, c in J['li5'].items()]
prods = [(nm.split('*'), F(c)) for nm, c in J['products'].items()]
assert all(len(n) == 5 for n, _ in prods if 'pi2' not in n)
T_coef, z5_coef = F(J['T_coef']), F(J['z5_coef'])


def chi17(a):
    a %= 17
    return 0 if a == 0 else (1 if pow(a, 8, 17) == 1 else -1)


def value(dps):
    mp.mp.dps = dps + 20
    s = mp.sqrt(17)
    lv = {'lE17': mp.log(4 + s), 'ls17': mp.log(s), 'lp': mp.log((3 + s) / 2), "lp'": mp.log(abs((3 - s) / 2)), 'pi2': mp.pi ** 2}
    L5 = sum(chi17(a) * mp.zeta(5, mp.mpf(a) / 17) for a in range(1, 17)) / mp.mpf(17) ** 5
    tot = T_coef.numerator * L5 / s / T_coef.denominator + z5_coef.numerator * mp.zeta(5) / z5_coef.denominator
    for (a, b), c in args:
        x = mp.mpf(a.numerator) / a.denominator + mp.mpf(b.numerator) / b.denominator * s
        tot += c.numerator * li5(x, dps + 10) / c.denominator
    for nms, c in prods:
        tot += c.numerator * mp.fprod(lv[n] for n in nms) / c.denominator
    return tot, L5, lv


r1, L5, _ = value(1350)
r2, _, _ = value(900)
print(f'(1) residual at 1350 digits: {mp.nstr(r1, 5)};  at 900 digits: {mp.nstr(r2, 5)};  L(chi17,5) = {mp.nstr(L5, 25)}')

# (2) factorisations by PSLQ at 120 digits
mp.mp.dps = 120
s = mp.sqrt(17)
lv = [mp.log(4 + s), mp.log(s), mp.log((3 + s) / 2), mp.log(abs((3 - s) / 2))]


def fac(x):
    v = mp.pslq([mp.log(abs(x))] + lv, maxcoeff=10**6, maxsteps=10**6)
    assert v and v[0] != 0, f'no factorisation for {x}'
    assert all(vi % v[0] == 0 for vi in v[1:]), f'non-integral factorisation {v}'
    return [-vi // v[0] for vi in v[1:]]


worst, vx, v1m = 0, [], []
for (a, b), c in args:
    x = mp.mpf(a.numerator) / a.denominator + mp.mpf(b.numerator) / b.denominator * s
    ex, e1 = fac(x), fac(1 - x)
    worst = max(worst, abs(mp.log(abs(x)) - mp.fsum(e * l for e, l in zip(ex, lv))),
                abs(mp.log(abs(1 - x)) - mp.fsum(e * l for e, l in zip(e1, lv))))
    vx.append(ex); v1m.append(e1)
print(f'(2) all {len(args)} arguments: x and 1-x are S-units over the 4 letters; worst |log - exponents.letters| = {mp.nstr(worst, 3)}')

# (3) exact symbol: Li5(x) -> -(1-x) (x) x (x) x (x) x x;  l1...l5 -> sum over orderings
sym = defaultdict(F)
for ((ab, c), ex, e1) in zip(args, vx, v1m):
    for t in itertools.product(range(4), repeat=5):
        w = e1[t[0]] * ex[t[1]] * ex[t[2]] * ex[t[3]] * ex[t[4]]
        if w:
            sym[t] -= c * w
idx = {n: k for k, n in enumerate(LET)}
for nms, c in prods:
    if 'pi2' in nms:
        continue  # pi^2 has zero symbol
    for perm in itertools.permutations(range(5)):
        sym[tuple(idx[nms[i]] for i in perm)] += c
nz = sum(1 for v in sym.values() if v)
print(f'(3) exact symbol: {nz} nonzero entries -> {"ZERO" if nz == 0 else "NONZERO"}')

# (4) stats
coefs = [c for _, c in args] + [c for _, c in prods] + [z5_coef]
print(f'(4) Li5 terms {len(args)}, products {len(prods)}, T_coef {T_coef} (17^5 = {17**5}), z5 {z5_coef}; '
      f'max |numerator| {max(abs(c.numerator) for c in coefs)}, lcm denominator {lcm(*[c.denominator for c in coefs])}')
