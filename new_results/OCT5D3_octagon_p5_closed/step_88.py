import os, sys, pickle, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
from oct5d3_basis import (atom_list, li32_symbol, B, H, R2, ONE, ZERO, Q, add, scale, LI5_ARGS)
from fractions import Fraction as F

t0 = time.time()
atoms = atom_list()
vis = [a for a in atoms if a[1]=='vis']
LI32 = [('Li32(b,1)',(B,ONE)),('Li32(1/sqrt2,1)',(H,ONE)),('Li32(-b,1)',(-B,ONE)),
        ('Li32(b,-1)',(B,-ONE)),('Li32(-b,-1)',(-B,-ONE)),('Li32(1/sqrt2,-1)',(H,-ONE)),
        ('Li32(1-1/sqrt2,1)',(ONE-H,ONE)),('Li32(2-sqrt2,1)',(2-R2,ONE)),('Li32(-b^2,1)',(-B*B,ONE))]
li32_syms = [li32_symbol(x,y) for _,(x,y) in LI32]
all_syms = [a[2] for a in vis] + li32_syms
labels_all = [a[0] for a in vis] + [n for n,_ in LI32]
words = sorted(set().union(*[set(s) for s in all_syms]))
M = [[s.get(w, ZERO) for s in all_syms] for w in words]
nr, nc = len(M), len(M[0])
A = [row[:] for row in M]
pivots = []; r = 0
for c in range(nc):
    piv = next((i for i in range(r, nr) if A[i][c] != ZERO), None)
    if piv is None: continue
    A[r], A[piv] = A[piv], A[r]
    pv = A[r][c]
    A[r] = [x/pv for x in A[r]]
    for i in range(nr):
        if i != r and A[i][c] != ZERO:
            f = A[i][c]
            A[i] = [A[i][j] - f*A[r][j] for j in range(nc)]
    pivots.append(c); r += 1
    if r == nr: break
print("rank:", r, "of", nr, "atoms:", nc, "kernel dim:", nc-r, "(%.1fs)" % (time.time()-t0))
free = [c for c in range(nc) if c not in pivots]
K61 = []
for f in free:
    v = [ZERO]*nc; v[f] = ONE
    for i, c in enumerate(pivots): v[c] = -A[i][f]
    K61.append(v)
# which kernel vectors involve the Li32's beyond K32?
for v in K61:
    li32_part = [(labels_all[52+j], v[52+j]) for j in range(9) if v[52+j] != ZERO]
    if li32_part:
        print("mixed kernel vector, Li32 part:", [(n, str(z)) for n, z in li32_part])
pickle.dump({'K61': K61, 'labels': labels_all, 'nclassical': 52}, open("kernel61.pkl","wb"))
print("saved kernel61.pkl")
