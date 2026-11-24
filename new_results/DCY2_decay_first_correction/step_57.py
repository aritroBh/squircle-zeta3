import sympy as sp
pi = sp.pi; I_ = sp.I
wb = (1-I_)/sp.sqrt(2); om = (1+I_)/sp.sqrt(2)
b = pi*wb/2

S1 = pi*(1+I_)/(6*sp.sqrt(2))
g1 = -3*wb/4
B2 = sp.Rational(1,2) + 3*pi*wb**2/4
A1 = 1 + pi*wb**2/2
gamma1 = A1**2/(4*b)

total = sp.simplify(S1 + g1 + B2/(2*b) + gamma1)
varpi1 = ((24+12*pi+pi**2) + I_*(24-12*pi+pi**2))/(48*sp.sqrt(pi))
target = sp.simplify(varpi1*sp.sqrt(2/pi))
print("S1+g1+B2/2b+gamma1 =", sp.expand(total))
print("varpi1*sqrt(2/pi)  =", sp.expand(target))
print("diff =", sp.simplify(total - target))
print("numeric:", sp.N(total, 30))

# also the moment-sum cancellation: -3T0/2 + 15/(16b) + 3T0/2 - 15 om/(4 pi) + 15 om/(8 pi) = 0 ?
T0 = A1/(2*b)
ms = sp.simplify(-3*T0/2 + 15/(16*b) + 3*T0/2 - 15*om/(4*pi) + 15*om/(8*pi))
print("moment-sum O(y) coeff =", ms)
