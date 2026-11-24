import pickle, subprocess
wd = '<project>/research/08_theory/TOW5R_review'
t5 = '<project>/research/08_theory/TOW5_tower_depth10'
# compare freshly-produced N5_proved.pkl with TOW5's
mine = pickle.load(open(wd+'/tow5_copy/N5_proved.pkl','rb'))
theirs = pickle.load(open(t5+'/N5_proved.pkl','rb'))
import sympy as sp
d, m = sp.symbols('d m')
nm = sp.sympify(mine['N5'], locals={'d':d,'m':m,'I':sp.I})
nt = sp.sympify(theirs['N5'], locals={'d':d,'m':m,'I':sp.I})
print('fresh N5 == TOW5 N5_proved.pkl:', sp.expand(nm-nt)==0)
rm = sp.sympify(mine['R5'], locals={'d':d,'m':m,'I':sp.I})
rt = sp.sympify(theirs['R5'], locals={'d':d,'m':m,'I':sp.I})
print('fresh R5 == theirs:', sp.cancel(rm-rt)==0)
# also R5_raw.pkl (the pre-simplification raw total?) vs my own total rebuild - check consistency
raw = pickle.load(open(t5+'/R5_raw.pkl','rb'))
print('R5_raw.pkl keys:', list(raw.keys()))
