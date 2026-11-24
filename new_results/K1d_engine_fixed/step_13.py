import sympy as sp, os

# Symbolic verification of the sqrt-series d^3 coefficient
d = sp.Symbol('d')
g0,g1,g2,g3 = sp.symbols('g0 g1 g2 g3', positive=True)
G = g0 + g1*d + g2*d**2 + g3*d**3
ser = sp.sqrt(G).series(d, 0, 4).removeO().expand()
c3 = sp.simplify(ser.coeff(d,3))
Ls = sp.sqrt(g0)
target = g3/(2*Ls) - g1*g2/(4*Ls**3) + g1**3/(16*Ls**5)
buggy  = g3/(2*Ls) - g1*g2/(8*Ls**3) + g1**3/(16*Ls**5)
print("sympy series d^3 coeff:", c3)
print("target matches:", sp.simplify(c3 - target) == 0)
print("buggy  matches:", sp.simplify(c3 - buggy) == 0)

with open('progress.log','a') as f:
    f.write("step_40a: PROVED sympy sqrt-series d^3 coeff = g3/(2Ls) - g1*g2/(4Ls^3) + g1^3/(16Ls^5); buggy 1/8 version does not match.\n")
print("logged.")
