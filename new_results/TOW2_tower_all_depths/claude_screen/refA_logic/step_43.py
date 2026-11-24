import sys, os
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
from towlib import *
from fractions import Fraction as Fr
from math import comb
import sympy as sp

def e_from_data(key, delta, s):
    i,j = key; dd = i+j
    base = (2*i+1) if i==j else dd
    return sum( u(i,j,base+t, delta).scale(Fr((-1)**(s-t)*comb(s,t))) for t in range(0,s+1) )

vals = []
for dd in range(4, 12):
    for key in keys_of_degree(dd):
        i,j = key
        if 4 > 2*dd-4: continue
        base = (2*i+1) if i==j else dd
        if base+1 > 14: continue
        e = e_from_data(key, 4, 1); e00 = e_from_data(key, 0, 0)
        R = e/e00
        vals.append((key, R))
for key, R in vals:
    print(key, R.p0, R.p2)
