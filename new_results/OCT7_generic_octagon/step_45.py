import sympy as sp
l2, l3, l5 = sp.symbols('l2 l3 l5')
ell1 = 2*l2 - l3
k1 = sp.Rational(1,3); k2 = sp.Rational(1,9)
expr = sp.log((k1-k2)/k1)
print("expr:", expr, repr(expr))
test = sp.expand(expr*ell1)
print("expanded:", test)
sub = test.subs(sp.log(sp.Rational(2,3)), l3-l2)
print("after subs:", sub)
# Now check engine's actual D5
import mixed_engine_gen as E
print("\nengine D5 =", E.D5)
