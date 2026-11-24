
import sys
import sympy as sp
u, v = sp.symbols('u v'); I = sp.I; S2 = sp.sqrt(2)
tau10 = (-1-I)/2
a = 1+S2; b = S2-1
DL1v = 1/(v-b) + 1/(v+a) - 1/v
g = sp.Rational(3,5) + I*sp.Rational(2,7)
G = g/(v-I)**3
prod = sp.together(G*DL1v)
ap = sp.apart(prod, v, extension=[S2, I])
print("apart terms:")
for t in sp.Add.make_args(sp.expand(ap)):
    print("   ", sp.radsimp(t))
# expected top: g*tau10/(v-I)^2
print("expected:", sp.radsimp(g*tau10), "/(v-I)^2")
