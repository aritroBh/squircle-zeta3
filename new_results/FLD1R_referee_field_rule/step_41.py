import sys, pickle
sys.path.insert(0,'.')
# (a) FLD1 k=3 Z-sum form == W1 template's mu_zsum.pkl (step_23's claim)
W1='<project>/research/10_frontier/W1_corner_family/B_taylor_pi4/'
tpl = pickle.load(open(W1+'mu_zsum.pkl','rb'))
mine = pickle.load(open('mu_zsum_K5.pkl','rb'))
print('template keys:', list(tpl.keys())[:6], '... mine k=3 keys:', list(mine[3].keys()))
import sympy as sp
n=sp.Symbol('n')
ok=True
for zkey, E in mine[3].items():
    zk3 = zkey[:2] if len(zkey)>=2 else zkey   # template likely uses (z2,z3) arity 2/3?
# inspect template structure
for k in list(tpl.keys())[:3]:
    print('template entry:', k, type(tpl[k]), str(tpl[k])[:200])
