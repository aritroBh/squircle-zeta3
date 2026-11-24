"""natural_extension_test.py -- NUMERICAL illustration (the general statement is the exact theorem in
function_symbol_check.py / b2b2_test.py).

Literal "natural extension" R_nat(c) of the right-hand side R of conj:W4: every constant of R is replaced by a
function of c that reduces to it at c0 = 1/sqrt2, using G1's standard arguments where possible:
   b -> c/(1+c),  -b -> (c-1)/c,  -b^2 -> -(1-c)/(1+c),  1/sqrt2 -> c,  1-1/sqrt2 -> 1-c,  2-sqrt2 -> 1/(1+c),
   1/2 -> 1/2,  L = log a -> log((1+c)/c),  lambda -> log 2,  sqrt2 -> 1/c  (variant B: sqrt2 -> 2c).
Compare with the weight-4 component W_4(c) of the corner coefficient M_3 (G1's B3, proved for all n here).
"""
import sys
import mpmath as mp
sys.path.insert(0, '.')
import hpl
from check_b3_split import R3, PSI, MPHI


def W4(c, N):
    s = mp.mpf(0)
    for terms, z, wt in ((R3, mp.mpf(1), 4), (PSI, mp.zeta(2), 2), (MPHI, mp.zeta(3), 1)):
        for e, w, v in terms:
            if len(w) == wt:
                s += z*mp.mpf(v.numerator)/v.denominator*c**e*hpl.value(w, c, N)
    return s


def li(n, x):
    return mp.re(mp.polylog(n, x))


def R_nat(c, r2):
    b = c/(1 + c); mb = (c - 1)/c; mb2 = -(1 - c)/(1 + c)
    L = mp.log((1 + c)/c); lam = mp.log(2); pi2 = mp.pi**2; z3 = mp.zeta(3)
    Li2b = li(2, b)
    t = ((1 - 573*r2)*li(4, b) - (1 + 573*r2)*li(4, mb) + 56*r2*li(4, mb2)
         + 17*r2*li(4, c) - (mp.mpf(37)/4 + 17*r2/16)*li(4, mp.mpf(1)/2)
         + (10 + 86*r2)*li(4, 1 - c) - (10 - 86*r2)*li(4, 1/(1 + c))
         + L*((40 + 7*r2)*li(3, b) - (40 - 7*r2)*li(3, mb) - 72*li(3, c) + (mp.mpf(63)/8 - 3*r2/16)*z3)
         + lam*((6 + 13*r2)*li(3, b) + (6 - 13*r2)*li(3, mb) - 47*r2/2*li(3, c) + (329*r2/128 - mp.mpf(3)/8)*z3)
         + 30*Li2b**2 + Li2b*((30 + 7*r2/2)*L**2 - 24*lam*L - 19*r2/8*lam**2 + (r2/24 - mp.mpf(19)/4)*pi2)
         + (mp.mpf(103)/384 - 4963*r2/5760)*pi2**2
         + pi2*((mp.mpf(67)/24 + 185*r2/96)*L**2 + (41*r2/12 - mp.mpf(35)/12)*lam*L + (mp.mpf(5)/12 - 1105*r2/384)*lam**2)
         + (mp.mpf(13)/24 - 25*r2/24)*L**4 - (mp.mpf(73)/6 + 29*r2/3)*L**3*lam + (mp.mpf(63)/8 + 67*r2/16)*L**2*lam**2
         + (mp.mpf(41)/24 + 17*r2/6)*L*lam**3 + (137*r2/192 - mp.mpf(61)/192)*lam**4)
    return t


if __name__ == '__main__':
    mp.mp.dps = 60
    c0 = 1/mp.sqrt(2)
    N = 520
    print('c          W4(c)                      R_nat(c), sqrt2->1/c - W4      R_nat(c), sqrt2->2c - W4')
    for c in [mp.mpf('0.05'), mp.mpf('0.3'), mp.mpf('0.5'), mp.mpf('0.65'), mp.mpf('0.7'), c0, mp.mpf('0.72'), mp.mpf('0.75')]:
        w = W4(c, N)
        dA = R_nat(c, 1/c) - w
        dB = R_nat(c, 2*c) - w
        print(f'{mp.nstr(c, 8):10s} {mp.nstr(w, 20):26s} {mp.nstr(dA, 12):30s} {mp.nstr(dB, 12)}')
    print('at c0 the literal extension equals R and W4(c0) equals the squircle W4 (check_b3_split.out).')
