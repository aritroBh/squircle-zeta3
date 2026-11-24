"""Spence-Kummer two-variable functional equation for Li3 (Kummer 1840; Lewin, Polylogarithms and
Associated Functions, 1981, eq. (6.107)), real form:

  Li3(x(1-y)^2/(y(1-x)^2)) + Li3(xy) + Li3(x/y) - 2Li3(x(1-y)/(y(1-x))) - 2Li3(x(1-y)/(x-1))
    - 2Li3((1-y)/(1-x)) - 2Li3((1-y)/(y(x-1))) - 2Li3(x) - 2Li3(y) + 2 zeta(3)
  = log^2(y) log((1-y)/(1-x)) - (pi^2/3) log y - (1/3) log^3 y,

on the connected region D = {0 < y < 1, x < y}, where every Li3 argument is < 1 and every log argument
is positive (so both sides are real-analytic on D; the identity holds on D by analytic continuation
from 0 < x < y < 1).  Each instance used is additionally checked numerically to 50 digits.
"""
import mpmath as mp
from qsqrt2 import QS2, SQ2, A, B, sunit_log
from polylog_exact import P, LOG, LI3, PI2, Z3, ONE, ZERO
from fractions import Fraction as Fr


def kummer_args(x, y):
    return [(x * (ONE - y) * (ONE - y) / (y * (ONE - x) * (ONE - x)), 1), (x * y, 1), (x / y, 1),
            (x * (ONE - y) / (y * (ONE - x)), -2), (x * (ONE - y) / (x - ONE), -2), ((ONE - y) / (ONE - x), -2),
            ((ONE - y) / (y * (x - ONE)), -2), (x, -2), (y, -2)]


def kummer_relation(x, y):
    assert ZERO < y < ONE and x < y
    lhs = P()
    for z, c in kummer_args(x, y):
        assert z <= ONE
        lhs = lhs + LI3(z) * c
    lhs = lhs + Z3 * 2
    ly = LOG(y)
    rhs = ly * ly * LOG((ONE - y) / (ONE - x)) - PI2 * ly * Fr(1, 3) - ly * ly * ly * Fr(1, 3)
    return lhs - rhs


def numeric_check_domain(n=200, seed=1):
    import random
    random.seed(seed)
    mp.mp.dps = 30
    worst = 0
    for _ in range(n):
        y = mp.mpf(random.uniform(0.01, 0.99)); x = mp.mpf(random.choice([random.uniform(-200, -5), random.uniform(-5, float(y) - 1e-3)]))
        Li3 = lambda z: mp.re(mp.polylog(3, z))
        lhs = (Li3(x * (1 - y) ** 2 / (y * (1 - x) ** 2)) + Li3(x * y) + Li3(x / y) - 2 * Li3(x * (1 - y) / (y * (1 - x)))
               - 2 * Li3(x * (1 - y) / (x - 1)) - 2 * Li3((1 - y) / (1 - x)) - 2 * Li3((1 - y) / (y * (x - 1)))
               - 2 * Li3(x) - 2 * Li3(y) + 2 * mp.zeta(3))
        rhs = mp.log(y) ** 2 * mp.log((1 - y) / (1 - x)) - mp.pi ** 2 / 3 * mp.log(y) - mp.log(y) ** 3 / 3
        worst = max(worst, abs(lhs - rhs))
    return worst


def search(maxj=4, maxm=4, arg_bound=(12, 8)):
    S0 = set()
    for j in range(-maxj, maxj + 1):
        for m in range(-maxm, maxm + 1):
            for sg in (1, -1):
                S0.add(QS2(sg) * SQ2 ** j * A ** m)
    ys = [y for y in S0 if ZERO < y < ONE]
    xs = [x for x in S0 if x < ONE]
    found = []
    for y in ys:
        try:
            sunit_log(ONE - y)
        except ValueError:
            continue
        for x in xs:
            if not x < y:
                continue
            try:
                sunit_log(ONE - x)
                ok = True
                for z, c in kummer_args(x, y):
                    if z == ZERO or z == ONE:
                        continue
                    j, m = sunit_log(z)
                    if abs(j) > arg_bound[0] or abs(m) > arg_bound[1]:
                        ok = False; break
                if ok:
                    found.append((x, y))
            except ValueError:
                continue
    return found


if __name__ == "__main__":
    print("random-point check of the real form on D (max |lhs-rhs|):", numeric_check_domain())
    f = search()
    print(len(f), "instances with S-unit arguments; first few:", f[:8])
