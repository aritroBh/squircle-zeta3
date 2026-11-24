import sys
sys.path.insert(0,'.')
from oct9c_common import *
D = load_all()
mp.mp.dps = 60

# 1) W3+Wlo == B_e ?
s123, s0, Wp0 = D['s123'], D['s0'], D['Wp0']
engs = {e: get_engine(e) for e in range(4)}
W3 = {1: s123[1]['W3'], 2: s123[2]['W3'], 3: s123[3]['W3'], 0: Wp0['W3_e0']}
Wlo = {1: s123[1]['Wlo'], 2: s123[2]['Wlo'], 3: s123[3]['Wlo'], 0: Wp0['Wlo_e0']}
for e in range(4):
    diff = sp.expand(D['brackets'][f'B{e}'] - (W3[e]+Wlo[e]))
    print(f"edge {e}: B-(W3+Wlo) = {diff}")

# 2) lengths: E_arc / B_e
P3 = D['P3']
for e in range(4):
    NUM = atom_values(engs[e], e, D['lam'], 60)
    Bval = mp.mpf(sp.N(D['brackets'][f'B{e}'].subs(list(NUM.items())), 60))
    E = mp.mpf(P3['vals'][e])
    print(f"edge {e}: E_arc/B = {mp.nstr(E/Bval, 25)}   sqrt(LEN2) = {mp.nstr(mp.sqrt(mp.mpf(str(LEN2[e]))), 25)}")
