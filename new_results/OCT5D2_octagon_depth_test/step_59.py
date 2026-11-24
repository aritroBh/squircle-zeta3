import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
from qsqrt2 import QS2 as Q
ZERO = Q(0)

data = pickle.load(open("symbol_fit2.pkl","rb"))
words, M, b, cols = data['words'], data['M'], data['b'], data['cols']
nr, nc = len(M), len(M[0])
A = [row[:] + [b[i]] for i, row in enumerate(M)]
pivots = []; r = 0
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
    pivots.append(c); r += 1
    if r == nr: break
assert r == 32
x = [ZERO]*nc
for i, c in enumerate(pivots): x[c] = A[i][nc]
from collections import defaultdict
byk = defaultdict(list)
maxnum = 0
for (label, kind), xv in zip(cols, x):
    if xv != ZERO:
        byk[kind].append((label, xv))
        maxnum = max(maxnum, max(abs(xv.p.numerator), abs(xv.q.numerator)))
for k in ['Li5','p41','p32','p23','log5','p32b','p221','li32']:
    print("== %s: %d nonzero ==" % (k, len(byk[k])))
    for label, xv in byk[k][:40]: print("   %-30s %s" % (label, xv))
print("max numerator magnitude:", maxnum, "(%d digits)" % len(str(maxnum)))
pickle.dump({'x': x}, open("symbol_fit2_solution.pkl","wb"))
