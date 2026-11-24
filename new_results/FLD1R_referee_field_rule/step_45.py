import pickle, json
from fractions import Fraction as Fr
store = pickle.load(open('hpl_constructed_K7.pkl','rb'))
# -1 usage in rational part; term stats
rat = store[(0,0,0,0,0,0)]
uses_minus1 = sum(1 for (e,w,q) in rat if -1 in w)
print(f'k=7 rational part: {len(rat)} terms; terms using letter -1: {uses_minus1} ({100*uses_minus1/len(rat):.1f}%)')
tot = sum(len(v) for v in store.values())
print('total k=7 HPL terms:', tot, 'over', len(store), 'zeta-monomials')
# mu decay sanity: mu_{7,n} * n^2 magnitude (rational part), n=10..30
MU = json.load(open('mu_exact_K7_N30.json'))
for nv in (10,20,30):
    v = Fr(MU['7'][str(nv)].get('0,0,0,0,0,0','0'))
    print(f'n={nv}: |mu_{{7,n}}[1]|*n^2 = {float(abs(v)*nv*nv):.4f}')
import os
print(sorted(f for f in os.listdir('.') if f.startswith('step_') and f.endswith('.py'))[-8:])
