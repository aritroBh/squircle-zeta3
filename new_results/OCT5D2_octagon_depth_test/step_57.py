import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
from qsqrt2 import QS2 as Q
ZERO = Q(0)

data = pickle.load(open("symbol_fit.pkl","rb"))
words, M, b, cols = data['words'], data['M'], data['b'], data['cols']
nr, nc = len(M), len(M[0])
A = [row[:] + [b[i]] for i, row in enumerate(M)]
pivots = []
r = 0
for c in range(nc):
    piv = None
    for i in range(r, nr):
        if A[i][c] != ZERO: piv = i; break
    if piv is None: continue
    A[r], A[piv] = A[piv], A[r]
    pv = A[r][c]
    A[r] = [x/pv for x in A[r]]
    for i in range(nr):
        if i != r and A[i][c] != ZERO:
            f = A[i][c]
            A[i] = [A[i][j] - f*A[r][j] for j in range(nc+1)]
    pivots.append(c)
    r += 1
    if r == nr: break
inconsistent_rows = [i for i in range(nr) if all(A[i][c]==ZERO for c in range(nc)) and A[i][nc]!=ZERO]
print("rank:", r, "of", nr, "rows;", "inconsistent rows:", len(inconsistent_rows))
print("S(W5_oct) in span of classical+Li32 symbols:", not inconsistent_rows)
if not inconsistent_rows:
    x = [ZERO]*nc
    for i, c in enumerate(pivots): x[c] = A[i][nc]
    from collections import defaultdict
    byk = defaultdict(list)
    for (label, kind), xv in zip(cols, x):
        if xv != ZERO: byk[kind].append((label, xv))
    for k in ['Li5','p41','p32','p23','log5','li32']:
        print("== %s: %d nonzero ==" % (k, len(byk[k])))
        for label, xv in byk[k]: print("   %-24s %s" % (label, xv))
    pickle.dump({'x': x, 'pivots': pivots}, open("symbol_fit_solution.pkl","wb"))
