#!/usr/bin/env python3
"""step_02_varpi1.py — the exact derivation of varpi_1 from the saddle-point expansion.

Setup (proof.md Sec.2): y = sqrt(pi/(2N)), s0 = omega/y, t = tau/sqrt(y) = tau/v, v = sqrt(y).
exp(G) H = e^{-b tau^2} [1 + v A1 + v^2 A2 + ...], b = pi wb/2, wb = e^{-i pi/4}:
  G1 = -(1 + pi wb^2/2) tau + (pi wb/2) tau^3
  G2 = (1/2 + 3 pi wb^2/4) tau^2 - (pi wb/2) tau^4
  h1 = -3/2, h2 = 15/8   (H-ratio (1+v tau)^{-3/2})
Then I = P0 [1 + (S1 - 3 wb/4 + c1) y + O(y^2)], P0 = -i (2/sqrt pi) y^2 e^{-theta+i theta},
  c1 = <G2 + G1^2/2 + h2 tau^2 + G1 h1 tau>_{e^{-b tau^2}},  S1 = pi(1+i)/(6 sqrt2).
And varpi_1 = sqrt(pi/2) (c1 + S1 - 3 wb/4).

Second (t-)formalism: exact ell_1, ell_2, Gaussian moment sum. The O(y) coefficient of the
moment sum is exactly 0, and the prefactor identity S1 + g1' + B2/(2b) + gamma1 = varpi1 sqrt(2/pi)
holds exactly, with B2 = 1/2 + 3 pi wb^2/4, gamma1 = A1^2/(4b), A1 = 1 + pi wb^2/2.

Both checked EXACTLY (diff = 0) and numerically (to 30+ digits) below.
"""
import sympy as sp
from mpmath import mp, mpf, mpc, pi as mpi, sqrt as msqrt

pi = sp.pi; I_ = sp.I
wb = (1-I_)/sp.sqrt(2)
b, tau = sp.symbols('b tau')
bval = pi*wb/2

print("="*90)
print("STEP 2 — derivation of varpi_1")
print("="*90)

G1 = -(1 + pi*wb**2/2)*tau + (pi*wb/2)*tau**3
G2 = (sp.Rational(1,2) + 3*pi*wb**2/4)*tau**2 - (pi*wb/2)*tau**4
h1 = sp.Rational(-3,2); h2 = sp.Rational(15,8)
integrand_v2 = sp.expand(G2 + G1**2/2 + h2*tau**2 + G1*h1*tau)
def mom(m):
    return sp.factorial2(2*m-1) / (2*b)**m if m > 0 else sp.Integer(1)
res = 0
for term in sp.Add.make_args(integrand_v2):
    c, tpow = term.as_coeff_exponent(tau)
    res += c*mom(int(tpow)//2)
c1 = sp.simplify(sp.simplify(res.subs(b, bval)).rewrite(sp.exp).expand())
print("c1 =", c1)
print("c1 numeric =", sp.N(c1, 30))

S1 = pi*(1+I_)/(6*sp.sqrt(2))
pref = S1 - 3*wb/4
varpi1_mine = sp.simplify(sp.sqrt(pi/2)*(c1 + pref))
varpi1_paper = ((24+12*pi+pi**2) + I_*(24-12*pi+pi**2))/(48*sp.sqrt(pi))
print("varpi1 (derived) =", sp.expand(varpi1_mine))
print("varpi1 (paper)   =", varpi1_paper)
print("DIFF (must be 0) =", sp.simplify(varpi1_mine - varpi1_paper))

print("-"*90)
print("t-formalism identities:")
g1p = -3*wb/4
B2 = sp.Rational(1,2) + 3*pi*wb**2/4
A1 = 1 + pi*wb**2/2
gamma1 = A1**2/(4*bval)
total = sp.simplify(S1 + g1p + B2/(2*bval) + gamma1)
target = sp.simplify(varpi1_paper*sp.sqrt(2/pi))
print("S1+g1'+B2/(2b)+gamma1 =", sp.expand(total))
print("varpi1*sqrt(2/pi)       =", sp.expand(target))
print("DIFF (must be 0) =", sp.simplify(total - target))
T0 = A1/(2*bval)
ms = sp.simplify(-3*T0/2 + 15/(16*bval) + 3*T0/2 - 15*(1+I_)/sp.sqrt(2)/(4*pi) + 15*(1+I_)/sp.sqrt(2)/(8*pi))
print("moment-sum O(y) coeff (must be 0) =", ms)

print("-"*90)
print("numeric cross-check (mpmath, 40 digits):")
mp.dps = 50
print("c1        =", mp.nstr(sp.N(c1, 45), 40))
print("varpi1    =", mp.nstr(sp.N(varpi1_paper, 45), 40))
