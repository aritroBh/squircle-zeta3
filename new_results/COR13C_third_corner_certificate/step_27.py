import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, '../HEXC1_hexagon_certificate')
sys.path.insert(0, '.')
import cmzv6
from cor13c_common import target4
from flint import nmod_mat, ctx
ctx.threads = 1
P = 1000000009
def cmod(c, p):
    return (c.numerator % p) * pow(c.denominator % p, -1, p) % p if hasattr(c,'numerator') else c % p
t0 = time.monotonic()
lower = {}
for n in (1,2,3):
    rows, _ = cmzv6.make_relations(n, lower)
    lower[n] = rows
cols4 = cmzv6.words(4)
rows4, counts4 = cmzv6.make_relations(4, lower)
print(f'weight 4: words={len(cols4)} rels={len(rows4)} families={counts4} ({time.monotonic()-t0:.1f}s)')
ids = {w:i for i,w in enumerate(cols4)}
mat = nmod_mat(len(rows4), len(cols4), P)
for i,row in enumerate(rows4):
    for w,c in row:
        j = ids.get(w)
        if j is not None: mat[i,j] = cmod(c,P)
print(f'matrix built ({time.monotonic()-t0:.1f}s), RREF...')
rr, r = mat.rref(inplace=True)
print(f'rank={r} quotient={len(cols4)-r} ({time.monotonic()-t0:.1f}s)')
T = target4()
assert set(T) <= set(cols4)
trow = [(w, c) for w, c in sorted(T.items())]
mat2 = nmod_mat(len(rows4)+1, len(cols4), P)
for i,row in enumerate(rows4+[trow]):
    for w,c in row:
        j = ids.get(w)
        if j is not None: mat2[i,j] = cmod(c,P)
r2 = mat2.rank()
print(f'target support={len(T)} rank+target={r2} member_modp={r2==r}')
# second prime for confidence
P2 = 1000000097
matB = nmod_mat(len(rows4), len(cols4), P2)
for i,row in enumerate(rows4):
    for w,c in row:
        j = ids.get(w)
        if j is not None: matB[i,j] = cmod(c,P2)
rB = matB.rank()
matB2 = nmod_mat(len(rows4)+1, len(cols4), P2)
for i,row in enumerate(rows4+[trow]):
    for w,c in row:
        j = ids.get(w)
        if j is not None: matB2[i,j] = cmod(c,P2)
rB2 = matB2.rank()
print(f'prime2: rank={rB} rank+target={rB2} member_modp2={rB2==rB}')
with open('progress.log','a') as f:
    f.write(f'step_27: weight-4 quotient={len(cols4)-r}, T4 member mod p = {r2==r} (two primes)\n')
