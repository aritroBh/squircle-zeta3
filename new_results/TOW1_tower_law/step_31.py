from towlib import *
import sympy as sp
from fastred import pf_to_sympy, pf_fold
v = sp.Symbol('v')
# test pf_fold on simple building blocks
def test_fold(P, name):
    Re = pf_to_sympy(P)
    ok = sp.simplify(sp.expand(pf_to_sympy(pf_fold(P)) - Re.subs(v,1/v)/v**2))==0
    print(name, ok)
    if not ok:
        print('  diff:', sp.simplify(sp.expand(pf_to_sympy(pf_fold(P)) - Re.subs(v,1/v)/v**2)))

Z, O, I2 = ZERO, ONE, IQ
test_fold([{BQ: O}, [], [], Z], 'simple pole b')
test_fold([{ZERO: O}, [], [], Z], 'simple pole 0')
test_fold([{}, [O], [], Z], 'order1 +i')
test_fold([{}, [Z, O], [], Z], 'order2 +i')
test_fold([{}, [Z, Z, O], [], Z], 'order3 +i')
test_fold([{}, [], [O], Z], 'order1 -i')
test_fold([{}, [], [Z, O], Z], 'order2 -i')
