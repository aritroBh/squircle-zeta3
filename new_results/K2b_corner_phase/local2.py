
import sympy as sp

# Local series of Btilde at z_+ = i alpha:
# Btilde = (E u - Ebar)/(1 + sqrt(2 i s E) sqrt(u) v), u = e^zeta - 1, v = sqrt(1 + E u/(2 i s))
E, Eb, s = sp.symbols('E Eb s')
zet = sp.Symbol('zet')
ORD = 8
u = sp.series(sp.exp(zet)-1, zet, 0, ORD+1).removeO()          # zet + zet^2/2 + ...
wser = u/zet                                                     # u/zet
sq_u = sp.sqrt(zet)*sp.series(sp.sqrt(wser), zet, 0, ORD).removeO()
v = sp.series(sp.sqrt(1 + E*u/(2*sp.I*s)), zet, 0, ORD).removeO()
a = sp.sqrt(2*sp.I*s*E)*sq_u/sp.sqrt(zet)*v    # a = a0*(1 + a1 zet + a2 zet^2 + ...) [series multiplier of sqrt(zet)]
a = sp.series(a, zet, 0, ORD).removeO()
# Denominator 1 + a sqrt(zet); expand 1/(1+X), X = a sqrt(zet)
X = a*sp.sqrt(zet)
inv = 0
sign = 1
Xp = 1
for j in range(0, 2*ORD+1):
    inv += sign*Xp
    Xp = sp.series(sp.expand(Xp*X), zet, 0, ORD+1).removeO()
    sign *= -1
inv = sp.series(sp.expand(inv), zet, 0, ORD+1).removeO()
Bser = sp.series(sp.expand((-Eb + E*u)*inv), zet, 0, ORD+1).removeO()
# extract coefficients of zet^{1/2}, zet^{3/2}, zet^{5/2}
sq = sp.sqrt(zet)
gam = sp.expand(Bser.coeff(sq,1))
g1 = sp.expand(Bser.coeff(sq,3)/gam)
g2 = sp.expand(Bser.coeff(sq,5)/gam)
b0 = sp.expand(Bser.coeff(sq,0))
b1 = sp.expand(Bser.coeff(sq,2))
print("b0 =", b0)
print("b1 =", b1)
print("gam =", gam, " (expect Eb*sqrt(2 i s E))")
print("g1 =", sp.simplify(g1))
print("g2 =", sp.simplify(g2))
# substitute Eb = 1/E and simplify with s free; also express in terms of E and s
g1s = sp.simplify(g1.subs(Eb, 1/E))
g2s = sp.simplify(g2.subs(Eb, 1/E))
print("g1(E,s) =", sp.expand(g1s))
print("g2(E,s) =", sp.expand(g2s))

g3 = sp.expand(Bser.coeff(sq,7)/gam)
g4 = sp.expand(Bser.coeff(sq,9)/gam)
print("g3 =", sp.simplify(g3.subs(Eb,1/E)))
print("g4 =", sp.simplify(g4.subs(Eb,1/E)))
