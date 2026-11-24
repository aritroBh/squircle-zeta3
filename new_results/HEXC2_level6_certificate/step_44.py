import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/HEXC1_hexagon_certificate')
sys.path.insert(0, '.')
import cmzv6
from hexc2_common import target2, target3
from flint import nmod_mat, ctx
ctx.threads = 1
P = 1000000009

def cmod(c, p):
    return (c.numerator % p) * pow(c.denominator % p, -1, p) % p if hasattr(c,'numerator') else c % p

def mat_modp(rows, cols, p=P):
    ids = {w:i for i,w in enumerate(cols)}
    mat = nmod_mat(len(rows), len(cols), p)
    for i,row in enumerate(rows):
        for w,c in row:
            j = ids.get(w)
            if j is not None: mat[i,j] = cmod(c,p)
    return mat

t0 = time.monotonic()
lower = {}
targets = {2: target2(), 3: target3()}
for n in (1,2,3):
    cols = cmzv6.words(n)
    rows, counts = cmzv6.make_relations(n, lower)
    mat = mat_modp(rows, cols)
    r = mat.rank()
    print(f"weight {n}: words={len(cols)} rels={len(rows)} rank={r} quotient={len(cols)-r} ({time.monotonic()-t0:.1f}s)")
    lower[n] = rows
    if n in targets:
        T = targets[n]
        trow = [(w, c) for w, c in sorted(T.items())]
        mat2 = mat_modp(rows + [trow], cols)
        r2 = mat2.rank()
        print(f"   target support={len(T)} rank+target={r2} member_modp={r2==r}")
