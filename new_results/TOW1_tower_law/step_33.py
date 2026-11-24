from towlib import *
import sympy as sp
from fastred import pf_to_sympy, pf_fold, pf_mul_DL, pf_scale, DL1, DL2, v as fv
def test_fold(P, name):
    Re = pf_to_sympy(P)
    d = sp.expand(pf_to_sympy(pf_fold(P)) - Re.subs(fv,1/fv)/fv**2)
    print(name, 'OK' if d==0 else f'FAIL {d}')
Z, O = ZERO, ONE
test_fold([{BQ: O}, [], [], Z], 'pole b')
test_fold([{ZERO: O}, [], [], Z], 'pole 0')
test_fold([{}, [O], [], Z], 'ord1 +i')
test_fold([{}, [Z, O], [], Z], 'ord2 +i')
test_fold([{}, [Z, Z, O], [], Z], 'ord3 +i')
test_fold([{}, [Z,Z,Z,Z,O], [], Z], 'ord5 +i')
test_fold([{BQ:O, ZERO:O}, [Z,O],[Z,Z,O], Z], 'mixed')
print('fold checks done', flush=True)
