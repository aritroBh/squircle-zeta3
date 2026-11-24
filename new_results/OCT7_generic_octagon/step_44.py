import sympy as sp
print("sp.log(Rational(2,3)) =", sp.log(sp.Rational(2,3)))
print("sp.log(Rational(5/2)) =", sp.log(sp.Rational(5,2)))
print("sp.log(Rational(2,5)) =", sp.log(sp.Rational(2,5)))
print(repr(sp.log(sp.Rational(2,3))))
# check what D5 etc contain
import mixed_engine_gen as E
for atom in sp.preorder_traversal(E.D5):
    if isinstance(atom, sp.log): print("D5 log atom:", atom)
for atom in sp.preorder_traversal(E.D7):
    if isinstance(atom, sp.log): print("D7 log atom:", atom)
for atom in sp.preorder_traversal(E.D8):
    if isinstance(atom, sp.log): print("D8 log atom:", atom)
print("sympy version:", sp.__version__)
