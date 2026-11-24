"""step_25.py -- PART 3a: M2 corner coefficient, t-integrand derivation (symbolic c).
M2 = int (1/g)[ A0(r2v - r1v log(1+g) + log^2(1+g)/2) + A0(r1v - log(1+g)) log g + (A0-1)log^2 g/2 ] dg
   = int [ RA2 log^2 g + RB2 log g log(1+g) + RC2 log^2(1+g) ] dt,  RA2,RB2,RC2 rational over Q(c).
"""
import sympy as sp
t, c = sp.symbols('t c')
ca = 2*c**2 - 1
c2a = sp.expand(2*ca**2 - 1)
r1, r2 = 2*c - 1, -1 - 2*c
P = t**2 - 2*ca*t + 1
g = (t**2 - 1)/(2*(t - ca))
J = P/((t**2 - 1)*(t - ca))            # (1/g) dg/dt
A0 = P/((t - r1)*(t - r2))
S0 = P**2/(4*(t - ca)**2)
# building blocks in t (proved identities from step_12)
cag_g2 = -(t**2 - 1)*(t**2 - 2*ca*t + c2a)/(4*(t - ca)**2)     # ca*g - g^2
g2_2_cag = (t**2 - 1)*(t**2 - ca*t + ca**2 - 1)/(2*(t - ca)**2) # 2g^2 - ca*g
R1c = sp.cancel(cag_g2/S0)                                     # (ca g-g^2)/S0
R2c = sp.cancel(g2_2_cag/(2*S0) - cag_g2**2/(2*S0**2))         # r2v/log^2 g
one_cag = sp.cancel((-ca*t**2 + 2*t - ca)/(2*(t - ca))/S0)     # (1-ca g)/S0... check
one_cag_direct = sp.cancel((1 - ca*g)/S0)
print("1-cag identity ok:", sp.simplify(one_cag - one_cag_direct) == 0)
print("R1c =", sp.factor(R1c))
print("R2c =", sp.factor(R2c))

RA2 = sp.cancel(J*(A0*(R2c + R1c + sp.Rational(1,2)) - sp.Rational(1,2)))
RB2 = sp.cancel(-J*A0*(1 + R1c))
RC2 = sp.cancel(J*A0*sp.Rational(1,2))
print("\nRA2 =", sp.factor(RA2))
print("\nRB2 =", sp.factor(RB2))
print("\nRC2 =", sp.factor(RC2))
import json
with open("step_25_out.json", "w") as f:
    json.dump({"RA2": str(sp.factor(RA2)), "RB2": str(sp.factor(RB2)), "RC2": str(sp.factor(RC2))}, f, indent=1)
print("\nsaved step_25_out.json")
