import sys, os, pickle
from math import gcd
from functools import reduce
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
import sympy as sp
from fractions import Fraction as Fr

with open(WORK+'/expansions.pickle','rb') as f: E = pickle.load(f)
freewords = [(7,b,c) for b in (4,5,6,7) for c in range(8)]

def row_lcm(coeffs):
    l = 1
    for q in coeffs:
        l = l*q.denominator//gcd(l, q.denominator)
    return l

rows_int = []
scales = []
for w in freewords:
    coeffs = E[w]['coeffs']
    l = row_lcm(coeffs)
    scales.append(l)
    rows_int.append([int(q*l) for q in coeffs])
print('row scales (lcm of denominators):')
for w,l in zip(freewords, scales): print(w, l)

# find pivot set of 27 independent rows mod p
p = 2**61 - 1
M = sp.Matrix(rows_int)
ind = M.echelon_form()  # over QQ might be slow; do mod p instead
Mm = M.applyfunc(lambda x: x % p)
# manual pivots via gf
Mrref, pivots = Mm.rref()  # sympy rref over QQ - careful; use domain
print('done')
