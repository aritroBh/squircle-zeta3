# Planner (W4 session) verification of ZQ5c's weight-3 Zagier formulas w3_D5.json and w3_D13.json:
# (1) value: T_coef*L(chi_D,3)/sqrt(D) + z3_coef*zeta(3) + sum c*Li3(x) + sum c*products = 0, with mpmath's polylog
#     (real part), own Hurwitz L-value and own letter logs, at two precisions;
# (2) S-unit factorisations of x and 1-x over the letters, found here by PSLQ (not read from ZQ5c);
# (3) exact symbol of the right-hand side in (Z^d)^{(x)3} is ZERO (Li3(x) -> -(1-x) (x) x (x) x; pi^2, zeta(3), T have
#     zero symbol; products -> all orderings).
import sys, json, itertools
from fractions import Fraction as F
from collections import defaultdict
import mpmath as mp


def parse(s, D):  # 'a+b*sqrtD' or '(a+b*sqrtD)' or '(b*sqrtD)' or '(2)'
    s = s.strip('()').replace('+-', '-').replace(' ', '')
    tag = f'*sqrt{D}'
    if tag not in s:
        return F(s), F(0)
    i = s.rfind(tag); head = s[:i]
    k = max(head.rfind('+'), head.rfind('-'))
    if k <= 0:
        return F(0), F(head)
    return F(head[:k]), F(head[k:])


def chi(D, a):
    a %= D
    if a == 0: return 0
    return 1 if pow(a, (D - 1) // 2, D) == 1 else -1  # D prime = 1 mod 4 here (5, 13)


def check(fn):
    J = json.load(open(fn)); D = J['D']
    names = J['names']; gens = [parse(g, D) for g in J['letters']]
    args = [(parse(x, D), F(c)) for x, c in J['li3']]  if isinstance(J['li3'], list) else [(parse(x, D), F(c)) for x, c in J['li3'].items()]
    prods = J['products'].items() if isinstance(J['products'], dict) else J['products']
    prods = [(nm.split('*'), F(c)) for nm, c in prods]
    T, z3 = F(J['T_coef']), F(J['z3_coef'])
    out = []
    for dps in (650, 850):
        mp.mp.dps = dps + 20
        s = mp.sqrt(D)
        val = lambda ab: mp.mpf(ab[0].numerator) / ab[0].denominator + mp.mpf(ab[1].numerator) / ab[1].denominator * s
        lv = {n: mp.log(abs(val(g))) for n, g in zip(names, gens)}; lv['pi2'] = mp.pi ** 2
        L3 = mp.fsum(chi(D, a) * mp.zeta(3, mp.mpf(a) / D) for a in range(1, D)) / mp.mpf(D) ** 3
        # ZQ5c's w3 JSONs use the convention  T_coef * L(chi_D,3)/sqrt(D) = z3*zeta(3) + sum Li3 + products
        # (the opposite side from d17_final.json); a first run with T on the same side gave residual 2*T*L/sqrt(D).
        tot = -T.numerator * L3 / s / T.denominator + z3.numerator * mp.zeta(3) / z3.denominator
        big = 0
        for ab, c in args:
            x = val(ab); big = max(big, abs(x))
            tot += c.numerator * mp.re(mp.polylog(3, x)) / c.denominator
        for nms, c in prods:
            tot += c.numerator * mp.fprod(lv[n] for n in nms) / c.denominator
        out.append((dps, tot))
    print(f'{fn}: D={D}, {len(args)} Li3, {len(prods)} products, T={T}, z3={z3}, max|x|={mp.nstr(big, 6)}')
    print('  (1) residuals: ' + ';  '.join(f'{d} dps: {mp.nstr(r, 4)}' for d, r in out) + f';  L(chi_{D},3)/sqrt{D} = {mp.nstr(L3 / mp.sqrt(D), 25)}')
    mp.mp.dps = 120
    s = mp.sqrt(D)
    val = lambda ab: mp.mpf(ab[0].numerator) / ab[0].denominator + mp.mpf(ab[1].numerator) / ab[1].denominator * s
    lg = [mp.log(abs(val(g))) for g in gens]

    def fac(x):
        v = mp.pslq([mp.log(abs(x))] + lg, maxcoeff=10 ** 6, maxsteps=10 ** 6)
        assert v and v[0] != 0 and all(vi % v[0] == 0 for vi in v[1:]), f'no integral factorisation for {x}: {v}'
        e = [-vi // v[0] for vi in v[1:]]
        assert abs(mp.log(abs(x)) - mp.fsum(a * l for a, l in zip(e, lg))) < mp.mpf(10) ** -100
        return e
    d = len(gens); sym = defaultdict(F)
    for ab, c in args:
        x = val(ab); ex, e1 = fac(x), fac(1 - x)
        for t in itertools.product(range(d), repeat=3):
            w = e1[t[0]] * ex[t[1]] * ex[t[2]]
            if w: sym[t] -= c * w
    idx = {n: k for k, n in enumerate(names)}
    for nms, c in prods:
        if 'pi2' in nms: continue
        for perm in itertools.permutations(range(3)):
            sym[tuple(idx[nms[i]] for i in perm)] += c
    nz = sum(1 for v in sym.values() if v)
    print(f'  (2) all {len(args)} arguments: x and 1-x are S-units over the {d} letters (PSLQ, checked to 1e-100)')
    print(f'  (3) exact symbol: {nz} nonzero entries -> {"ZERO" if nz == 0 else "NONZERO"}')


for fn in sys.argv[1:]:
    check(fn)
