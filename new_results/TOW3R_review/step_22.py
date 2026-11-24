import os, pickle, sys, datetime
from fractions import Fraction as Fr
base = '<project>/research/08_theory/'
work = base + 'TOW3R_review/'
os.makedirs(work, exist_ok=True)
os.chdir(work)
H = pickle.load(open(base+'TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
keys = [k for k in H if k != '_meta']
print(len(keys), 'keys')
print('meta:', H.get('_meta'))
# inspect one key
K = (3,1)
print('deltas for (3,1):', sorted(H[K].keys()))
print('levels for (3,1) delta 0:', sorted(H[K][0].keys())[:10], '...', sorted(H[K][0].keys())[-3:])
print('sample:', H[K][0][4], H[K][5][4] if 5 in H[K] and 4 in H[K][5] else None)
