import pickle
from fractions import Fraction as Fr
base = '<project>/research/08_theory/'
H = pickle.load(open(base+'TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
K=(3,2)
print('deltas (3,2):', sorted(H[K].keys()))
print('levels delta 5:', sorted(H[K][5].keys()))
print('u(3,2,5,0)=',H[K][0][5])
print('u(3,2,5,5)=',H[K][5][5])
print('u(3,2,5,6)=',H[K][6][5])
# diagonal key (2,2) degree 4 -> deltas?
print('deltas (2,2):', sorted(H[(2,2)].keys()))
K=(3,3)
print('deltas (3,3):', sorted(H[K].keys()), 'levels d6:', sorted(H[K][6].keys())[:6])
print('u(3,3,6,0)=',H[K][0][6])
print('u(3,3,6,6)=',H[K][6][6])
