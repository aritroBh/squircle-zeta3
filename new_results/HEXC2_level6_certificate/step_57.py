import sys, time, pickle, json
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/HEXC1_hexagon_certificate')
sys.path.insert(0, '.')
import cmzv6
from fractions import Fraction as Fr
import exact_linalg as XL

# (a) is zeta(2) = pi^2/6 in the weight-2 span?  vector: (-1,0) + 3*(1,1)
lower = {}
for n in (1,2):
    rows, _ = cmzv6.make_relations(n, lower)
    lower[n] = rows
cols2 = cmzv6.words(2)
v = {(-1,0): Fr(1), (1,1): Fr(3)}
ids = {w:i for i,w in enumerate(cols2)}
from flint import nmod_mat
P = 1000000009
def rank_of(rows):
    M = XL.rows_to_nmod(rows, len(cols2), ids, P)
    return M.rank()
r0 = rank_of(rows)
r1 = rank_of(rows + [[(w,c) for w,c in sorted(v.items())]])
print('zeta(2) = pi^2/6 in weight-2 span:', r0 == r1)

# (b) y heights
y4 = pickle.load(open('w4_y.pkl','rb'))
import math
h = [max(len(str(abs(c.numerator))), len(str(c.denominator))) for c in y4]
nz = sum(1 for c in y4 if c != 0)
print(f'w4 multipliers: {len(y4)} total, {nz} nonzero, max height {max(h)} digits, median {sorted(h)[len(h)//2]}')
cert3 = json.load(open('cert_w23.json'))
for n in ('2','3'):
    ys = [Fr(a)/Fr(b) for a,b in cert3[n]['y']]
    ys = [Fr(int(a), int(b)) for a,b in cert3[n]['y']]
    hn = [max(len(str(abs(c.numerator))), len(str(c.denominator))) for c in ys]
    print(f'w{n} multipliers: {len(ys)}, max height {max(hn)} digits')

# (c) dump weight-4 certificate json
pv = pickle.load(open('w4_pivots.pkl','rb'))
out = dict(cols=[list(w) for w in cmzv6.words(4)], prow=pv['prow'],
           y=[str(c) for c in y4])
json.dump(out, open('cert_w4.json','w'))
import os
print('cert_w4.json', os.path.getsize('cert_w4.json'), 'bytes')
