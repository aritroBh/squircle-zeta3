import pickle, sympy as sp, sys
work = '<project>/research/08_theory/TOW5_tower_depth10'
sys.path.insert(0, work)
d, m = sp.symbols('d m')
Lraw = pickle.load(open(work+'/L_simplified.pkl','rb'))
print("keys:", sorted(Lraw.keys()))
L8 = sp.sympify(Lraw[8], locals={'d':d,'m':m,'I':sp.I})
print("L(8) cached:")
print(L8)
nn, dd = sp.fraction(sp.together(L8))
print("L8 den factored:", sp.factor(dd))
# check the TOW3R sqrt2-aware step_04
import os
p04 = '<project>/research/08_theory/TOW3R_review/step_04.py'
print("="*80)
print(open(p04).read())
