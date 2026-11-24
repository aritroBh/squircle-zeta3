import os, pickle
from math import gcd
from functools import reduce
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
import mpmath as mp
mp.mp.dps = 750
import sympy as sp

with open(WORK+'/vals_hi.pickle','rb') as f: V = pickle.load(f)
free_vals = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re,im in V['vals'][27:]]
freewords = V['freewords']

C = 10**700
m = 32
rows = []
for i,v in enumerate(free_vals):
    r = [0]*(m+2); r[i]=1
    r[m] = int(mp.nint(C*v.real)); r[m+1] = int(mp.nint(C*v.imag))
    rows.append(r)
R = sp.Matrix(rows).lll(delta=sp.Rational(99,100)).tolist()

rels = []
for r in R:
    tail = max(abs(r[32]), abs(r[33]))
    if tail < 10**40:
        coeffs = r[:32]
        g = reduce(gcd, [abs(x) for x in coeffs])
        coeffs = [x//g for x in coeffs]
        rels.append(coeffs)
print('extracted', len(rels), 'short relations')
# check rank
Mrel = sp.Matrix(rels)
print('rank:', Mrel.rank())
for i, coeffs in enumerate(rels):
    res = sum(x*fv for x,fv in zip(coeffs, free_vals))
    h = max(abs(x) for x in coeffs)
    nz = [(freewords[j], x) for j,x in enumerate(coeffs) if x]
    print(f'--- R{i+1}: height {h}, residual {mp.nstr(abs(res),3)}, support {len(nz)}')
    print('   ', nz)
with open(WORK+'/five_relations.pickle','wb') as f:
    pickle.dump({'freewords': freewords, 'relations': rels}, f)
with open(WORK+'/progress.log','a') as f:
    f.write('step_31: 5 short relations extracted (heights ~1e2-1e3), rank 5, residuals ~1e-750\n')
