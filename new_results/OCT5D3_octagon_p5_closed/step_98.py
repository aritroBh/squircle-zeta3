import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
from fractions import Fraction as F

ak = pickle.load(open("atom_kernel.pkl","rb"))
rels, names = ak['rels'], ak['names']
print("atom names (%d):" % len(names), names[:6], "...")
# relation vectors are 102-dim (p,q per atom). Keep only genuinely short/exact ones (residuals were ~0 for the first ones)
vecs = [v for h, v in rels if h < 10**6]
print("using %d short relations" % len(vecs))

# preference order: K combos first (eliminate), then MIX, then K32, then planner
def pref(nm):
    if nm.startswith("K") and nm != "K32": return 0
    if nm.startswith("MIX"): return 1
    if nm == "K32": return 2
    return 3
order = sorted(range(len(names)), key=lambda i: pref(names[i]))
colperm = []
for i in order:
    colperm += [2*i, 2*i+1]
inv = [0]*102
for newc, oldc in enumerate(colperm): inv[oldc] = newc

M = [[F(row[c]) for c in colperm] for row in vecs]
nr, nc = len(M), 102
A = [row[:] for row in M]
pivots = []; r = 0
for c in range(nc):
    piv = next((i for i in range(r, nr) if A[i][c] != 0), None)
    if piv is None: continue
    A[r], A[piv] = A[piv], A[r]
    pv = A[r][c]
    A[r] = [x/pv for x in A[r]]
    for i in range(nr):
        if i != r and A[i][c] != 0:
            f = A[i][c]
            A[i] = [A[i][j]-f*A[r][j] for j in range(nc)]
    pivots.append(c); r += 1
print("independent relations:", r)
pivot_atoms = sorted(set(colperm[c]//2 for c in pivots))
print("atoms eliminated (%d):" % len(pivot_atoms), [names[i] for i in pivot_atoms])
keep = [i for i in range(len(names)) if i not in pivot_atoms]
print("basis atoms kept (%d):" % len(keep), [names[i] for i in keep])
pickle.dump({'keep': keep, 'names': names}, open("basis_keep.pkl","wb"))
print("saved basis_keep.pkl")
