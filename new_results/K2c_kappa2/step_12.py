
import sympy as sp

# beta2 as exact polynomial in q (q = sqrt(Z)), g1, g2  [from step_11]
q, g1, g2 = sp.symbols('q g1 g2')
beta2 = (-g1*q**4/8 + 21*g1*q**2/8 + 9*g1/4 + 15*g2*q**2/4
         + q**6/288 - 5*q**4/48 + 65*q**2/96 + sp.Rational(15,8))

# local series coefficients (K2b local.py/local2.py, PROVED)
a  = sp.Symbol('a', positive=True)      # alpha
E  = sp.exp(sp.I*a)
s  = sp.sin(a)
g1f = -sp.Rational(3,4) - sp.I*E/(4*s)
g2f = sp.Rational(29,96) + E**2*(sp.Rational(1,2) + 1/(32*s**2)) + (sp.I*E/s)*(E**2/4 - sp.Rational(3,16))

Z = sp.I*a
qf = sp.sqrt(Z)                          # principal sqrt of i*a
b2 = beta2.subs({q: qf, g1: g1f, g2: g2f})
b2 = sp.expand(b2)

# kappa2(alpha) = conjugate(beta2(i alpha)); alpha real > 0
k2 = sp.conjugate(b2)
k2r = sp.simplify(sp.expand(sp.re(k2)))
k2i = sp.simplify(sp.expand(sp.im(k2)))
print("Re kappa2 =", k2r)
print("Im kappa2 =", k2i)
print()
# rewrite in terms of cot / csc for presentation
cot = sp.cot(a); csc = 1/sp.sin(a)
k2r2 = sp.simplify(sp.trigsimp(sp.expand(k2r.rewrite(sp.sin).rewrite(sp.cos))))
k2i2 = sp.simplify(sp.trigsimp(sp.expand(k2i.rewrite(sp.sin).rewrite(sp.cos))))
print("Re kappa2 (trigsimp) =", k2r2)
print("Im kappa2 (trigsimp) =", k2i2)
print()
# numeric values at the 5 angles
import mpmath
mpmath.mp.dps = 50
for name, av in [('pi/3', mpmath.pi/3), ('pi/2', mpmath.pi/2), ('2pi/3', 2*mpmath.pi/3),
                 ('3pi/4', 3*mpmath.pi/4), ('5pi/6', 5*mpmath.pi/6)]:
    f = sp.lambdify(a, k2, 'mpmath')
    val = f(av)
    print('kappa2(%s) = %s' % (name, mpmath.nstr(val, 30)))
print()
print('pi/2 closed pieces: 3/4 + pi^2/96 =', mpmath.nstr(mpmath.mpf(3)/4 + mpmath.pi**2/96, 30))
print('                    7*pi/96 - pi^3/2304 =', mpmath.nstr(7*mpmath.pi/96 - mpmath.pi**3/2304, 30))
