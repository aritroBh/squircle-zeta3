import sys, os, pickle
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
import mpmath as mp
mp.mp.dps = 750
import sympy as sp
from fractions import Fraction as Fr

with open(WORK+'/vals_hi.pickle','rb') as f: V = pickle.load(f)
basis_vals = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re,im in V['vals'][:27]]
free_vals = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re,im in V['vals'][27:]]
freewords = V['freewords']
with open(WORK+'/expansions.pickle','rb') as f: E = pickle.load(f)

# 32 x 27 rational matrix
M = sp.Matrix([[E[w]['coeffs'][j] for j in range(27)] for w in freewords])
print('rank:', M.rank())  # expect 27
ns = M.T.nullspace()  # left nullspace: vectors x (len 32) with x^T M = 0 -> wait, M is 32x27; left null = null(M^T)
print('left nullity:', len(ns))  # expect 5

def mag(x):
    from math import gcd
    from functools import reduce
    den = 1
    for q in x:
        den = den*q.denominator // __import__('math').gcd(den, q.denominator)
    ints = [int(q*den) for q in x]
    g = reduce(__import__('math').gcd, ints)
    ints = [a//g for a in ints]
    return ints, max(abs(a) for a in ints)

rels = []
for i, v in enumerate(ns):
    ints, mx = mag(list(v))
    rels.append(ints)
    print(f'R{i+1}: max|c|={mx}')
    print('   ', ints)
# verify residuals numerically
for i, ints in enumerate(rels):
    res = sum(c*fv for c,fv in zip(ints, free_vals))
    print(f'R{i+1} residual: {mp.nstr(abs(res),3)}  (height {max(abs(a) for a in ints)})')
with open(WORK+'/five_relations.pickle','wb') as f:
    pickle.dump({'freewords': freewords, 'relations': rels}, f)
with open(WORK+'/progress.log','a') as f:
    f.write('step_24: extracted 5 relations among the 32 free words (left nullspace of 32x27 expansion matrix); residuals printed\n')
