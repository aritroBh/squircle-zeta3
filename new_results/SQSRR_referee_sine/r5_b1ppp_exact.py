# Part 3: exact assembly of b_1'''(2) from Parseval at third order,
#   b_1''' = N'''/2 - 3 b_1' b_1'' - 3 sum_{n>=3} b_n' b_n'' ,  sum_{n>=3} b_n'b_n'' = log2*Sbeta2 + pi^2*SbetaC + R,
# and of b_1'' from Parseval at second order, b_1'' = N''/2 - b_1'^2 - Sbeta2.  N from the Gamma formula.
import sympy as sp, pickle
from dshuffle import ZS, l2
ex = pickle.load(open('r3_R_exact.pkl', 'rb'))
p, e = sp.symbols('p e')
N = 2*sp.sin(sp.pi/p)/sp.pi*sp.gamma(3/p)*sp.gamma(1-1/p)/sp.gamma(1+2/p)
Ns = sp.series(N.subs(p, 2+e), e, 0, 4).removeO()
toS = {sp.log(2): l2, sp.pi: sp.sqrt(6*ZS[2]), sp.zeta(3): ZS[3]}
def clean(x):
    x = sp.expand(sp.expand_func(sp.expand_log(x, force=True)))
    x = x.subs({sp.polygamma(1, sp.Rational(3, 2)): sp.pi**2/2-4, sp.polygamma(2, sp.Rational(3, 2)): -14*sp.zeta(3)+16})
    x = sp.expand(sp.expand_log(x, force=True))
    return sp.expand(x.subs(toS).subs(ZS[2]**2, sp.Rational(5, 2)*ZS[4]))
N2 = clean(2*Ns.coeff(e, 2)); N3 = clean(6*Ns.coeff(e, 3))
print("N''(2)  =", N2); print("N'''(2) =", N3)
b1p = (l2-1)/2
b1pp = sp.expand(N2/2 - b1p**2 - ex['sum_{n>=3} beta^2'])
print("b_1''(2) =", b1pp)
T3 = sp.expand(l2*ex['sum_{n>=3} beta^2'] + 6*ZS[2]*ex['sum beta C'] + ex['R'])
b1ppp = sp.expand(N3/2 - 3*b1p*b1pp - 3*T3)
b1ppp = sp.expand(b1ppp.subs(ZS[2]**2, sp.Rational(5, 2)*ZS[4]))
print("b_1'''(2) =", b1ppp)
pi = sp.sqrt(6*ZS[2])
C4 = (-sp.Rational(3, 8) + sp.Rational(33, 16)*l2 - 3*pi**2/16 - sp.Rational(9, 8)*l2**2 + sp.Rational(33, 64)*ZS[3]
      - pi**2*l2/64 + l2**3/8 + 3*pi**4/512 - sp.Rational(93, 256)*ZS[5] + sp.Rational(9, 256)*pi**2*ZS[3] - pi**4*l2/512)
C4 = sp.expand(sp.expand(C4).subs(ZS[2]**2, sp.Rational(5, 2)*ZS[4]))
print("b_1''' - C4 =", sp.simplify(b1ppp - C4))
num = {l2: sp.log(2), **{ZS[s]: sp.zeta(s) for s in range(2, 8)}}
print("b_1'''(2) =", sp.N(b1ppp.subs(num), 70))
# odd-weight bookkeeping, in the basis pi^2 zeta3, zeta3, zeta5 (pi^2 = 6 zeta2)
for nm, part in (("N'''/2", N3/2), ("-3 pi^2 Sum beta C", -3*6*ZS[2]*ex['sum beta C']), ("-3R", -3*ex['R']), ('total', b1ppp)):
    part = sp.expand(part)
    print('  %-20s zeta3: %-8s  pi^2 zeta3: %-8s  zeta5: %s' % (nm, part.coeff(ZS[3]).subs({ZS[2]: 0, l2: 0}),
          part.coeff(ZS[3]).coeff(ZS[2])/6, part.coeff(ZS[5])))
