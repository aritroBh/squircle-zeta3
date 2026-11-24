import sys, pickle, os
from fractions import Fraction as Fr
TH = '<project>/research/08_theory'
sys.path.insert(0, TH + '/WTS_weight_k_structure/theory')
def kpath(k):
    if k <= 18: return f'{TH}/WTS_weight_k_structure/theory/kern_{k}.pkl'
    return f'{TH}/TWOEND1_weight_poles/red/kern_{k}.pkl'
_c = {}
def kern(k):
    if k not in _c:
        _c[k] = pickle.load(open(kpath(k), 'rb'))
    return _c[k]
KEYS5 = [(5,0),(4,1),(3,2)]
POLES = {'0':(0,0),'a':(1,1),'-a':(-1,-1),'-b':(1,-1)}
def comp(k, key, pole):
    """4-vector of Fractions of the real residue at pole (r0,r1); zeros if absent"""
    r = kern(k)[key]['real']
    return tuple(r.get(pole, (Fr(0),)*4))
def poleset(k, key):
    return set(kern(k)[key]['real'].keys())
