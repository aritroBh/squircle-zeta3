import sympy as sp
# T_c involution check
t, c = sp.symbols('t c')
T = c*(1-t)/(c+(1-c)*t)
TT = sp.simplify(T.subs(t, T))
print("T_c(T_c(t)) =", TT, " residual:", sp.simplify(TT - t))

# c values at the Landen pairs for m=3,4,5:
# m=4 octagon: x1=1-sqrt(2), x2 = x1/(x1-1); c = x2/(1+x1)
s2 = sp.sqrt(2); x1_4 = 1-s2; x2_4 = sp.simplify(x1_4/(x1_4-1)); c4 = sp.simplify(x2_4/(1+x1_4))
print("octagon: x1 =", x1_4, " x2 =", x2_4, " c =", c4)
# m=5 decagon: x1 = -1/phi
phi = (1+sp.sqrt(5))/2; x1_5 = -1/phi; x2_5 = sp.simplify(x1_5/(x1_5-1)); c5 = sp.simplify(x2_5/(1+x1_5))
print("decagon: x1 =", x1_5, " x2 =", x2_5, " c =", c5)
# m=3 hexagon: guess x1 = -1/2
x1_3 = sp.Rational(-1,2); x2_3 = sp.simplify(x1_3/(x1_3-1)); c3 = sp.simplify(x2_3/(1+x1_3))
print("hexagon candidate: x1 =", x1_3, " x2 =", x2_3, " c =", c3)
