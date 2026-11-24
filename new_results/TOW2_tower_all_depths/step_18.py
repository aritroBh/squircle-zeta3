import sys, os
sys.dont_write_bytecode = True
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
from towlib import *
from math import comb

def e_from_data(key, delta, s):
    """Newton coefficient e(delta,s): Delta^s u at k=d (off-diag); for diag use ehat via k=2p+1 base."""
    i,j = key; d = i+j
    if i == j:
        p = i
        return sum( u(i,j,2*p+1+t, delta).scale(Fr((-1)**(s-t)*comb(s,t))) for t in range(0,s+1) )
    return sum( u(i,j,d+t, delta).scale(Fr((-1)**(s-t)*comb(s,t))) for t in range(0,s+1) )

# print the full e-table for several keys, normalized by e(0,0)
for key in [(5,0),(4,1),(3,2),(3,3),(6,0)]:
    i,j = key; d = i+j
    e00 = e_from_data(key,0,0)
    print(f'--- key {key}, d={d}, e00={e00}')
    maxd = 2*d-4
    for delta in range(0, maxd+1):
        row = []
        for s in range(0, delta//2+1):
            # check data availability: need k = d+s <= 14 (off) or 2p+1+s<=14
            kk = (d+s) if i>j else (2*i+1+s)
            if kk > 14: 
                row.append('n/a'); continue
            e = e_from_data(key, delta, s)
            r = e/e00 if not e.is_zero() else ZERO
            row.append(f's={s}: {r}')
        print(f'  d{delta}: ' + ' | '.join(row))
