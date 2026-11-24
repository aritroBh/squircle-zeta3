import pickle, sympy as sp
work = '<project>/research/08_theory/TOW5_tower_depth10'
d, m = sp.symbols('d m')
p5 = pickle.load(open(work+'/N5_proved.pkl','rb'))
N5 = sp.sympify(p5['N5'], locals={'d':d,'m':m,'I':sp.I})
print("N_5 =")
print(sp.expand(N5))
# also sl_5 closed form denominator
import sys; sys.path.insert(0, work)
import tow5lib as T5
sn, sd = sp.fraction(sp.together(T5.sl5))
print("sl_5 den factored:", sp.factor(sd))
print("sl_5 num degree in d:", sp.degree(sp.expand(sn), d))
