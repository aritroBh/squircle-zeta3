"""step_12.py -- GENERAL-ANGLE THEOREM, part A: symbolic derivation with c = cos(a/2) free.
[1] substitution identities (Euler t = g + sqrt(S0)) proved symbolically;
[2] t-integrand = Astar(t)*log g + Bstar(t)*log(1+g) proved symbolically;
[3] universal residues: simplify(res - claimed) == 0 for all 5 poles;
[4] full partial-fraction identities simplify to 0.
"""
import sympy as sp

t, c = sp.symbols('t c')
ca = 2*c**2 - 1              # cos a
c2a = sp.expand(2*ca**2 - 1) # cos 2a = 8c^4 - 8c^2 + 1
r1, r2 = 2*c - 1, -1 - 2*c
P = t**2 - 2*ca*t + 1
D = (t**2 - 1)*(t - ca)*(t - r1)*(t - r2)

print("[1] substitution identities (g = (t^2-1)/(2(t-ca)))")
g = (t**2 - 1)/(2*(t - ca))
S0 = 1 + g**2 - 2*ca*g
dgdt = sp.diff(g, t)
ids = {
 "1+g = (t-r1)(t-r2)/(2(t-ca))": sp.simplify(1 + g - (t - r1)*(t - r2)/(2*(t - ca))),
 "S0  = P^2/(4(t-ca)^2)":       sp.simplify(S0 - P**2/(4*(t - ca)**2)),
 "dg/dt = P/(2(t-ca)^2)":       sp.simplify(dgdt - P/(2*(t - ca)**2)),
 "1-ca*g = (-ca t^2+2t-ca)/(2(t-ca))": sp.simplify(1 - ca*g - (-ca*t**2 + 2*t - ca)/(2*(t - ca))),
 "ca*g-g^2 = -(t^2-1)(t^2-2ca t+c2a)/(4(t-ca)^2)":
     sp.simplify(ca*g - g**2 + (t**2 - 1)*(t**2 - 2*ca*t + c2a)/(4*(t - ca)**2)),
 "cos2a = 8c^4-8c^2+1": sp.simplify(c2a - (8*c**4 - 8*c**2 + 1)),
}
for k, v in ids.items():
    print(f"   {k:55s} -> {v}")

print()
print("[2] t-integrand: M1 = int (1/g)[ A0*((1-ca g)/S0) log g - A0 log(1+g) - log g ] dg,")
print("    A0 = sqrt(S0)/(1+g),  sqrt(S0) = +P/(2(t-ca)) on path (both factors >0).")
lg, l1g = sp.symbols('lg l1g')   # formal log g, log(1+g)
A0 = (P/(2*(t - ca)))/((t - r1)*(t - r2)/(2*(t - ca)))
A0 = sp.cancel(A0)
S0v = P**2/(4*(t - ca)**2)
Ig = (1/g)*( A0*((1 - ca*g)/S0v)*lg - A0*l1g - lg )
integ = sp.cancel(Ig*dgdt)
X = sp.cancel(integ.coeff(lg)); Y = sp.cancel(integ.coeff(l1g))
Astar = P**2/D - P/((t**2 - 1)*(t - ca)) - (t**2 - 2*ca*t + c2a)/((t - ca)*(t - r1)*(t - r2))
Bstar = -P**2/D
print("   simplify(coeff of log g     - Astar) =", sp.simplify(X - Astar))
print("   simplify(coeff of log(1+g)  - Bstar) =", sp.simplify(Y - Bstar))

print()
print("[3] universal residues (symbolic c)")
poles = [sp.Integer(1), sp.Integer(-1), ca, r1, r2]
names = ["1", "-1", "cos a", "-1+2c", "-1-2c"]
expA = [sp.Integer(0), sp.Integer(-2), sp.Integer(1), -c, c]
expB = [sp.Integer(-1), sp.Integer(1), sp.Integer(-1), 2*c, -2*c]
def prod_over(rho, xs):
    p = sp.Integer(1)
    for x in xs:
        if x != rho: p *= (rho - x)
    return p
ok = True
for rho, nm, ea, eb in zip(poles, names, expA, expB):
    resA = P.subs(t, rho)**2/prod_over(rho, poles)
    if rho in (1, -1, ca):
        resA -= P.subs(t, rho)/prod_over(rho, [x for x in (sp.Integer(1), sp.Integer(-1), ca)])
    if rho in (ca, r1, r2):
        resA -= (rho**2 - 2*ca*rho + c2a)/prod_over(rho, [x for x in (ca, r1, r2)])
    resB = -P.subs(t, rho)**2/prod_over(rho, poles)
    dA = sp.simplify(resA - ea); dB = sp.simplify(resB - eb)
    ok &= (dA == 0 and dB == 0)
    print(f"   rho={nm:>6}: resA-({ea}) = {dA},   resB-({eb}) = {dB}")
print("   all residues universal:", ok)

print()
print("[4] full partial fractions")
pfA = sum(ea/(t - rho) for rho, ea in zip(poles, expA))
pfB = sum(eb/(t - rho) for rho, eb in zip(poles, expB))
print("   simplify(Astar - sum a_rho/(t-rho)) =", sp.simplify(Astar - pfA))
print("   simplify(Bstar - sum b_rho/(t-rho)) =", sp.simplify(Bstar - pfB))
