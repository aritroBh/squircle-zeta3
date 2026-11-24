import pickle, sys
import sympy as sp
d='../OCT9b_octagon_debug/'
def L(name):
    return pickle.load(open(d+name,'rb'))
brackets = L('brackets.pkl')
for k,v in brackets.items():
    print('====', k, type(v))
    if isinstance(v, dict):
        for kk in v: print('   ', kk, type(v[kk]))
    else:
        print(v)
