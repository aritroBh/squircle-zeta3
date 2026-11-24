# (1) Exact proof of sum_{n>=3} b_n' C_n = (1 - zeta3)/64 and sum_{n>=3} b_n'^2 (rational-function sums, Hurwitz reduction).
# (2) Check of the identified closed forms of R and b_1''' at 250 digits against both Euler-Maclaurin runs.
import sympy as sp
from symsum import i, sum_inf, ZETA, LN2, GAM
from mpmath import mp, mpf, log, pi, zeta, nstr
k = i   # summation index k >= 1, n = 4k + eps, a = 4k, c = 4k + 2 eps
tot_C = 0; tot_sq = 0
for eps in (1, -1):
    n = 4*k + eps; a = 4*k; c = 4*k + 2*eps
    beta = sp.Rational(2)*eps/(a*c**2)
    C = -eps*n/(4*a**2*c)                     # pi^2 coefficient of b_n'' (= -2 eps n/(a^2 c) * 1/8)
    tot_C += sum_inf(sp.together(beta*C), 1)
    tot_sq += sum_inf(sp.together(beta**2), 1)
print('sum b_n\' C_n  =', sp.nsimplify(sp.expand(tot_C)))
print('sum b_n\'^2    =', sp.expand(tot_sq), '  (zeta2 = pi^2/6, zeta4 = pi^4/90)')
mp.dps = 260
z2, z3, z5, l2 = zeta(2), zeta(3), zeta(5), log(2)
Rcf = mpf(5)/16 - 13*pi**2/96 + 47*z3/64 + pi**2*l2/32 - pi**4/384 + 31*z5/256 + pi**2*z3/256
Bcf = (-mpf(3)/8 + 33*l2/16 - 3*pi**2/16 - 9*l2**2/8 + 33*z3/64 - pi**2*l2/64 + l2**3/8 + 3*pi**4/512
       - 93*z5/256 + 9*pi**2*z3/256 - pi**4*l2/512)
for fn in ('R_260_300.txt', 'R_260_400.txt'):
    v = dict(l.split() for l in open(fn))
    print(fn, ' R - closed form =', nstr(mpf(v['R'])-Rcf, 3), "  b1''' - closed form =", nstr(mpf(v['b1ppp'])-Bcf, 3))
print('R  =', nstr(Rcf, 60)); print("b1''' =", nstr(Bcf, 60))
