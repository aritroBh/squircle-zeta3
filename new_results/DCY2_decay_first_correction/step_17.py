import sympy as sp
pi = sp.pi; I_ = sp.I
wb = (1-I_)/sp.sqrt(2)

# c1 from previous step (integral correction)
c1 = (1+I_)*(15-pi**2)/(8*sp.sqrt(2)*pi)
# prefactor: S1 = pi(1+i)/(6 sqrt2), and H series: -3 wb/4
S1 = pi*(1+I_)/(6*sp.sqrt(2))
pref = S1 - 3*wb/4

varpi1_mine = sp.simplify(sp.sqrt(pi/2)*(c1 + pref))
varpi1_paper = ((24+12*pi+pi**2) + I_*(24-12*pi+pi**2))/(48*sp.sqrt(pi))
print("mine  =", sp.simplify(varpi1_mine))
print("paper =", varpi1_paper)
print("diff  =", sp.simplify(varpi1_mine - varpi1_paper))
print("numeric mine :", sp.N(varpi1_mine, 40))
print("numeric paper:", sp.N(varpi1_paper, 40))
