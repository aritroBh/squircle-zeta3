import sys, os, pickle
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
import mpmath as mp
mp.mp.dps = 620
from fractions import Fraction as Fr

with open(WORK+'/basis_vals.pickle','rb') as f: B = pickle.load(f)
with open(WORK+'/free_vals_hi.pickle','rb') as f: F = pickle.load(f)
basis_vals = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re,im in B['vals']]
free_vals = {w: mp.mpc(mp.mpf(re), mp.mpf(im)) for w,(re,im) in F['freevals'].items()}
with open(WORK+'/expansions.pickle','rb') as f: E = pickle.load(f)

w = (7,4,0)
e = E[w]
# relation: a0*w + sum a_j b_j = 0  where coeffs[j] = -a_j/a0  => reconstruct integer vector
a0 = e['a0']
expr = e['coeffs']  # Fr list: w = sum expr[j] b_j
# value check: w - sum expr[j] b_j
res = free_vals[w] - sum(c*b for c,b in zip(expr, basis_vals))
print('w - sum expr*b  residual:', mp.nstr(abs(res),5))
print('height of expr denominators:', max(c.denominator for c in expr))
print('max numerator:', max(abs(c.numerator) for c in expr))
# also check the raw integer relation: a0 w + sum a_j b_j with a_j = -expr[j]*a0 (integers)
res2 = a0*free_vals[w] + sum((-c*a0)*b for c,b in zip(expr, basis_vals))
print('integer relation residual:', mp.nstr(abs(res2),5))
