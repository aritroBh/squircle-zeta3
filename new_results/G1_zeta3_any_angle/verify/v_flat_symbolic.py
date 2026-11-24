"""v_flat_symbolic.py -- verifier's own derivation of the zeta content of mu_{k,n} (Theorem A).

Integrand factor (y-form, checked by hand in verify.md):  E = exp(u*Ea) * (1 + y(e^{u D} - 1))^m,
  Ea = -(n-1) a + n b,  D = b - a,  m = 1 - 2n,  a = log y, b = log(1-y).
Expanded by explicit binomial series (no sympy.series of a log), coefficients of u^2, u^3.
mu_{k,n} = -(1/(n(2n-1))) sum_{j,i,l} e_{j,i,l} (n)_j/(2n)_j D_{i,l}(n+j, n),  D_{i,l} = d_p^i d_q^l B / B.
zeta parts: D_{1,1} -> -z2 ;  D_{2,1}, D_{1,2} -> 2 z3 + 2 kappa_{11} kappa_{10 or 01} -> 2 z3 - 2 z2 F, where
  F = psi(p) - psi(p+q) (for D_{2,1}) or psi(q) - psi(p+q) (for D_{1,2}).
Outputs: symbolic-n z2-coefficient of mu_2, z3-coefficient of mu_3; exact (n <= 40) z2-coefficient of mu_3
compared with (i) the solver's table and (ii) the Taylor coefficients of the solver's Psi(c) formula.
"""
import json
from fractions import Fraction as Fr
import sympy as sp

u, y, a, b, n, m = sp.symbols('u y a b n m')
Ea = -(n - 1)*a + n*b
D = b - a
X = y*(u*D + u**2*D**2/2 + u**3*D**3/6)
binom = 1 + m*X + m*(m - 1)/2*X**2 + m*(m - 1)*(m - 2)/6*X**3
expE = 1 + u*Ea + u**2*Ea**2/2 + u**3*Ea**3/6
E = sp.expand(binom*expE)
E2 = sp.Poly(sp.expand(E.coeff(u, 2).subs(m, 1 - 2*n)), y, a, b)
E3 = sp.Poly(sp.expand(E.coeff(u, 3).subs(m, 1 - 2*n)), y, a, b)


def poch(x, j):
    return sp.Mul(*[x + t for t in range(j)]) if j > 0 else sp.Integer(1)


# zeta(2) coefficient of mu_{2,n}, symbolic n
tot = 0
for (j, i, l), co in E2.terms():
    if (i, l) == (1, 1):
        tot += co*poch(n, j)/poch(2*n, j)*(-1)
z2c = sp.factor(sp.simplify(-tot/(n*(2*n - 1))))
print('zeta(2)-coefficient of mu_{2,n} (symbolic n):', z2c)
print('   degrees (i+l) of monomials in [u^2]:', sorted({i + l for (j, i, l) in E2.monoms()}))
# zeta(3) coefficient of mu_{3,n}, symbolic n
tot = 0
for (j, i, l), co in E3.terms():
    if (i, l) in ((2, 1), (1, 2)):
        tot += co*poch(n, j)/poch(2*n, j)*2
z3c = sp.factor(sp.simplify(-tot/(n*(2*n - 1))))
print('zeta(3)-coefficient of mu_{3,n} (symbolic n):', z3c)
print('   degrees (i+l) of monomials in [u^3]:', sorted({i + l for (j, i, l) in E3.monoms()}))


# zeta(2) coefficient of mu_{3,n}, exact for n <= 40
def H(N):
    return sum((Fr(1, t) for t in range(1, N + 1)), Fr(0))


terms3 = [(j, i, l, co) for (j, i, l), co in E3.terms()]
z2mu3 = {}
for nv in range(1, 41):
    tot = Fr(0)
    for (j, i, l, co) in terms3:
        if (i, l) not in ((2, 1), (1, 2)):
            continue
        cv = sp.Rational(co.subs(n, nv))
        cv = Fr(int(cv.p), int(cv.q))
        pr = Fr(1)
        for t in range(j):
            pr *= Fr(nv + t, 2*nv + t)
        p, q = nv + j, nv
        F = (H(p - 1) - H(p + q - 1)) if (i, l) == (2, 1) else (H(q - 1) - H(p + q - 1))
        tot += cv*pr*(-2*F)
    z2mu3[nv] = -tot/(nv*(2*nv - 1))
tab = json.load(open('flat_exact_N40_K3.json'))
ok_tab = all(Fr(tab['3'][str(nv)].get('(1, 0, 0)', '0')) == z2mu3[nv] for nv in range(1, 41))
print('zeta(2)-coefficient of mu_{3,n}: verifier == solver exact table for n<=40:', ok_tab)
print('   first values:', [str(z2mu3[k]) for k in range(1, 5)])


# Taylor coefficients of the solver's Psi(c) (result.md B3), exact Laurent series in c
NT = 2*40 + 6


def ser_H(word):
    a_ = [Fr(0)]*(NT + 3)
    last = word[-1]
    for k in range(1, NT + 3):
        a_[k] = Fr(1, k) if last == 1 else Fr((-1)**(k + 1), k)
    for letter in reversed(word[:-1]):
        b_ = [Fr(0)]*(NT + 3)
        if letter == 0:
            for k in range(1, NT + 3):
                b_[k] = a_[k]/k
        else:
            acc = Fr(0)
            for mm in range(1, NT + 2):
                acc = (-acc + a_[mm]) if letter == -1 else (acc + a_[mm])
                b_[mm + 1] = acc/(mm + 1)
        a_ = b_
    return a_


def shift(s, e):          # multiply series by c^e (e may be negative); returns dict power->coef
    return {k + e: v for k, v in enumerate(s) if v != 0}


def add(dst, src, f):
    for k, v in src.items():
        dst[k] = dst.get(k, Fr(0)) + f*v


H1, Hm1 = ser_H((1,)), ser_H((-1,))
W = [ser_H(w) for w in ((0, 1), (0, -1), (1, -1), (-1, 1))]
comb = [W[0][k] + W[1][k] + W[2][k] - W[3][k] for k in range(NT + 3)]
Psi = {0: Fr(1, 2)}
add(Psi, shift(H1, 0), Fr(3, 4)); add(Psi, shift(H1, -1), Fr(-1, 4)); add(Psi, shift(H1, -2), Fr(-1, 2))
add(Psi, shift(Hm1, 0), Fr(-3, 4)); add(Psi, shift(Hm1, -1), Fr(-1, 4)); add(Psi, shift(Hm1, -2), Fr(1, 2))
add(Psi, shift(comb, -1), Fr(1, 4)); add(Psi, shift(comb, 1), Fr(-1, 4))
bad_neg = {k: v for k, v in Psi.items() if k <= 0 and v != 0}
odd = {k: v for k, v in Psi.items() if 0 < k <= 2*40 and k % 2 == 1 and v != 0}
okPsi = all(Psi.get(2*nv, Fr(0)) == z2mu3[nv] for nv in range(1, 41))
print('Psi(c) series: nonzero coefficients at c^k, k<=0:', bad_neg, '; odd powers (k<=80):', len(odd))
print('Taylor coefficients of solver Psi(c) == verifier zeta(2)-coeff of mu_{3,n} for n<=40:', okPsi)
