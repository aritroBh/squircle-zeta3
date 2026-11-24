"""step_27.py -- PART 3c: exact partial fractions of RA2, RB2, RC2 over Q(c).
RA2 denom = 2(t+1)(t-r1)(t-r2)(t-ca) P^2  -> ansatz with simple real poles + (f1 t+f2)/P + (f3 t+f4)/P^2.
RB2 denom = (t-1)(t+1)(t-r1)(t-r2); RC2 denom = 2(t-1)(t+1)(t-r1)(t-r2)(t-ca): simple real poles only.
Solve linear systems over Q(c) exactly.
"""
import sympy as sp
t, c = sp.symbols('t c')
ca = 2*c**2 - 1
c2a = sp.expand(2*ca**2 - 1)
r1, r2 = 2*c - 1, -1 - 2*c
P = t**2 - 2*ca*t + 1
g = (t**2 - 1)/(2*(t - ca))
J = P/((t**2 - 1)*(t - ca))
A0 = P/((t - r1)*(t - r2))
S0 = P**2/(4*(t - ca)**2)
cag_g2 = -(t**2 - 1)*(t**2 - 2*ca*t + c2a)/(4*(t - ca)**2)
g2_2_cag = (t**2 - 1)*(t**2 - ca*t + ca**2 - 1)/(2*(t - ca)**2)
R1c = sp.cancel(cag_g2/S0)
R2c = sp.cancel(g2_2_cag/(2*S0) - cag_g2**2/(2*S0**2))
RA2 = sp.cancel(J*(A0*(R2c + R1c + sp.Rational(1,2)) - sp.Rational(1,2)))
RB2 = sp.cancel(-J*A0*(1 + R1c))
RC2 = sp.cancel(J*A0*sp.Rational(1,2))

def pf_simple(R, poles_list, denom):
    """residues at simple poles: N(rho)/prod_{sig != rho}(rho-sig) using known factorization."""
    num = sp.fraction(sp.cancel(R))[0]
    res = {}
    for rho, fac in poles_list:
        others = sp.cancel(denom/fac)
        res[rho] = sp.simplify(num.subs(t, rho)/others.subs(t, rho))
    return res

# RB2: denom (t-1)(t+1)(t-r1)(t-r2)
DB = (t - 1)*(t + 1)*(t - r1)*(t - r2)
polesB = [(sp.Integer(1), t - 1), (sp.Integer(-1), t + 1), (r1, t - r1), (r2, t - r2)]
resB = pf_simple(RB2, polesB, DB)
print("RB2 residues at (1,-1,r1,r2):")
for rho, _ in polesB:
    print(f"   rho={str(rho):>10}: {sp.factor(resB[rho])}")
chkB = sp.simplify(RB2 - sum(resB[rho]/(t - rho) for rho, _ in polesB))
print("   RB2 - sum res/(t-rho) =", chkB)

# RC2: denom 2(t-1)(t+1)(t-r1)(t-r2)(t-ca)
DC = 2*(t - 1)*(t + 1)*(t - r1)*(t - r2)*(t - ca)
polesC = polesB + [(ca, t - ca)]
resC = pf_simple(RC2, polesC, DC)
print("\nRC2 residues at (1,-1,r1,r2,ca):")
for rho, _ in polesC:
    print(f"   rho={str(rho):>10}: {sp.factor(resC[rho])}")
chkC = sp.simplify(RC2 - sum(resC[rho]/(t - rho) for rho, _ in polesC))
print("   RC2 - sum res/(t-rho) =", chkC)

# RA2: ansatz
DA = 2*(t + 1)*(t - r1)*(t - r2)*(t - ca)*P**2
NA = sp.fraction(sp.cancel(RA2))[0]
e1, e2, e3, e4, f1, f2, f3, f4 = sp.symbols('e1 e2 e3 e4 f1 f2 f3 f4')
ans = (e1/(t + 1) + e2/(t - r1) + e3/(t - r2) + e4/(t - ca)
       + (f1*t + f2)/P + (f3*t + f4)/P**2)
expr = sp.expand(ans*DA - NA)
coeffs = sp.Poly(expr, t).all_coeffs()
sol = sp.solve(coeffs, [e1, e2, e3, e4, f1, f2, f3, f4], dict=True)
print("\nRA2 partial fractions solution:")
for k, v in sol[0].items():
    print(f"   {k} = {sp.factor(v)}")
chkA = sp.simplify(RA2 - ans.subs(sol[0]))
print("   RA2 - partial fractions =", chkA)
