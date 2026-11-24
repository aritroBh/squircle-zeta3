import sys, time, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/HEXC1_hexagon_certificate')
sys.path.insert(0, '.')
import cmzv6
from hexc2_common import target4
from flint import nmod_mat, ctx
ctx.threads = 1
P = 1000000009
def cmod(c, p):
    return (c.numerator % p) * pow(c.denominator % p, -1, p) % p if hasattr(c,'numerator') else c % p

t0 = time.monotonic()
# rebuild lower weights (fast) then weight 4
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
trow = [(w, c) for w, c in sorted(T.items())]
mat2 = nmod_mat(len(rows4)+1, len(cols4), P)
for i,row in enumerate(rows4+[trow]):
    for w,c in row:
        j = ids.get(w)
        if j is not None: mat2[i,j] = cmod(c,P)
r2 = mat2.rank()
print(f'target support={len(T)} rank+target={r2} member_modp={r2==r}')
with open('progress.log','a') as f:
    f.write(f'step_44: quotient dims 3,8,23,{len(cols4)-r}; T2,T3 member mod p; T4 member mod p = {r2==r} ({time.monotonic()-t0:.0f}s)\n')
