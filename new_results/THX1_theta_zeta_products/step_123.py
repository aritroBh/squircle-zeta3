
import sys, time, pickle, itertools
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, flint
import numpy as np

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

rels_data = pickle.load(open('rels9_dsdual.pkl','rb'))
rels = [{k: Fr(n, d) for k, (n, d) in r.items()} for fam, tag, r in rels_data]
cols9 = pickle.load(open('cols9.pkl','rb'))
st = pickle.load(open('struct_w9.pkl','rb'))
B = st['indep_rows']; colpos = st['colpos']; pivcols = st['pivcols']
p0 = 4611686018427387847
r = len(B)
posmap = {c: jj for jj, c in enumerate(pivcols)}
free_cols = [c for c in range(len(cols9)) if c not in posmap]
freepos = {c: jj for jj, c in enumerate(free_cols)}
L = pickle.load(open('w9_L.pkl','rb'))['L']
arr = np.load('w9_Ainv.npy')
Ainv = flint.nmod_mat(r, r, p0)
for i in range(r):
    for j in range(r):
        Ainv[i,j] = int(arr[i,j])
RBf = flint.nmod_mat(r, 55, p0)
for i, bi in enumerate(B):
    for k, v in rels[bi].items():
        j = freepos.get(colpos[k])
        if j is not None:
            RBf[i, j] = (v.numerator % p0) * pow(v.denominator % p0, -1, p0) % p0

def free_reduce(combo):
    b = flint.nmod_mat(r, 1, p0)
    for k, v in combo.items():
        j = posmap.get(colpos[k])
        if j is not None:
            b[j, 0] = (v.numerator % p0) * pow(v.denominator % p0, -1, p0) % p0
    for i in range(r):
        b[i, 0] = (int(b[i,0]) * (L[i] % p0)) % p0
    x = Ainv * b
    xtRBf = x.transpose() * RBf
    out = [int(xtRBf[0, c]) for c in range(55)]
    for k, v in combo.items():
        j = freepos.get(colpos[k])
        if j is not None:
            out[j] = (out[j] - (v.numerator % p0) * pow(v.denominator % p0, -1, p0) % p0) % p0
    return [(-x) % p0 for x in out]

fstar = freepos[colpos[(0,1,1,1,1,1,1,0,-1)]]
print("X9f free index:", fstar)
# search {0,-1} words with nonzero fstar component
hits = []
for bits in itertools.product((0,-1), repeat=7):
    w = bits + (0,-1) if False else None
    # first letter in {0,-1}, last letter -1
for first in [0,-1]:
    for mid in itertools.product((0,-1), repeat=7):
        w = (first,) + mid + (-1,)
        red = free_reduce({w: Fr(1)})
        if red[fstar] % p0 != 0:
            hits.append((w, red[fstar] % p0))
print("{0,-1} words with nonzero X9f-component:", len(hits))
for w, v in hits[:12]:
    print(w, v)
pickle.dump({'hits': [w for w,_ in hits]}, open('w9_x9f_reps.pkl','wb'))
log("searched {0,-1} reps of X9f class: %d hits" % len(hits))
