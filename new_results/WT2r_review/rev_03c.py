
import sys
import sympy as sp
u, v = sp.symbols('u v'); I = sp.I; S2 = sp.sqrt(2)
tau10 = (-1-I)/2; tau20 = (-1+I)/2
a = 1+S2; b = S2-1
DL1v = 1/(v-b) + 1/(v+a) - 1/v
DL2v = 1/(v-a) + 1/(v+b) - 1/v

def top_coeff_at_i(expr, order):
    """exact coefficient of (v-I)^(-order) in rational expr via apart"""
    ap = sp.apart(sp.together(expr), v, extension=[S2, I])
    tot = 0
    for t in sp.Add.make_args(sp.expand(ap)):
        tt = sp.together(t)
        num, den = sp.fraction(tt)
        # den should be c*(v-I)^order
        if sp.simplify(den - (v-I)**order)==0:
            tot += num
    return sp.radsimp(sp.simplify(tot))

# Laurent coefficients of DL1, DL2 at v=I: c_n
import math
for nm, DL in [('DL1',DL1v),('DL2',DL2v)]:
    c1 = sp.simplify(sp.diff(DL,v).subs(v,I))
    c2 = sp.simplify(sp.diff(DL,v,2).subs(v,I)/2)
    print(f"{nm} at +i: c1={c1}, c2={sp.radsimp(c2)}   (c1 = tau: {sp.simplify(c1-tau10)==0 or sp.simplify(c1-tau20)==0})")

# B1 exact: random-ish rational instances of G
from fractions import Fraction as Fr
ok = True
for dd in range(2,13):
    M = 2*dd-1
    for trial in range(3):
        gv = [sp.Rational(2*j+1+trial, j+2) + I*sp.Rational(j-1+trial, j+3) for j in range(1, M+1)]
        G = sum(gv[j-1]/(v-I)**j for j in range(1, M+1))
        prod = G*DL1v
        tc = top_coeff_at_i(prod, 2*dd-2)
        ok &= (sp.simplify(sp.radsimp(tc - tau10*gv[M-1])) == 0)
print("(B1) [u^{-(2d-2)}](G*DL1) = tau1_0*[top of G], d=2..12, 3 exact instances each:", ok)
# same for DL2
ok2 = True
for dd in range(2,13):
    M = 2*dd-1
    gv = [sp.Rational(j+1, j+2) + I*sp.Rational(j, j+3) for j in range(1, M+1)]
    G = sum(gv[j-1]/(v-I)**j for j in range(1, M+1))
    tc = top_coeff_at_i(G*DL2v, 2*dd-2)
    ok2 &= (sp.simplify(sp.radsimp(tc - tau20*gv[M-1])) == 0)
print("(B1') [u^{-(2d-2)}](G*DL2) = tau2_0*[top of G], d=2..12:", ok2)
# and that LOWER coefficients of G do NOT contribute (set top=0 -> result 0):
ok3 = True
for dd in range(2,9):
    M = 2*dd-1
    gv = [sp.Rational(j+1, j+2) + I*sp.Rational(j, j+3) for j in range(1, M)] + [0]
    G = sum(gv[j-1]/(v-I)**j for j in range(1, M+1))
    tc = top_coeff_at_i(G*DL1v, 2*dd-2)
    ok3 &= (sp.simplify(tc) == 0)
print("(B1'') lower terms of G contribute nothing to top of product:", ok3)
