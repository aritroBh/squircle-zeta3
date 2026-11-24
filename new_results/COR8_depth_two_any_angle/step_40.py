
import pickle, numpy as np
from fractions import Fraction as Fr
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve81.pkl','rb'))
li4v, prodv, Rv, bwords = D['li4v'], D['prodv'], D['Rv'], D['bwords']

# ---- exact fraction Gaussian elimination helpers ----
def rref(rows, ncols):
    """rows: list of lists of Fr. returns (basis rows, pivots) in reduced row echelon form"""
    M = [list(r) for r in rows]
    basis = []
    pivots = []
    for r in M:
        for b, pw in zip(basis, pivots):
            if r[pw]:
                f = r[pw]
                r = [x - f*y for x, y in zip(r, b)]
        nz = next((j for j, x in enumerate(r) if x), None)
        if nz is None: continue
        f = r[nz]
        r = [x/f for x in r]
        # back-eliminate
        for i, (b, pw) in enumerate(zip(basis, pivots)):
            if b[nz]:
                g = b[nz]
                basis[i] = [x - g*y for x, y in zip(b, r)]
        basis.append(r); pivots.append(nz)
    return basis, pivots

# ---- product basis (choose 63 nice products) ----
names = list(prodv.keys())
vecs = list(prodv.values())
p = (1<<31)-1
def niceness(nm):
    simple = ['c)', '-c)', '1 - c)', '1/(c + 1)', '(1 - c)/(c + 1)', 'c/(c + 1)', '(c - 1)/(c + 1)',
              'c**2', '2*c/(c + 1)', '(1 - c)*(c + 1)']
    score = 0
    for i, s in enumerate(simple):
        if f"({s}" in nm or f"*{s}" in nm:
            score += 10 - i
    score -= nm.count('*') * 2 + len(nm)/100
    if nm.startswith('log'): score += 5
    return -score
order = sorted(range(len(names)), key=lambda k: niceness(names[k]))
basis_idx = []
B_rows = []  # mod-p rref rows
pivots = []
for k in order:
    v = np.array([x % p for x in vecs[k]], dtype=np.int64)
    for br, pw in zip(B_rows, pivots):
        if v[pw]:
            v = (v - int(v[pw])*br) % p
    nz = np.nonzero(v)[0]
    if len(nz):
        pw = int(nz[0])
        B_rows.append((v*pow(int(v[pw]), p-2, p)) % p)
        pivots.append(pw); basis_idx.append(k)
print("product basis size:", len(basis_idx))
Bn = [names[k] for k in basis_idx]
# exact 63x63 pivot matrix and its inverse
Bvecs = [vecs[k] for k in basis_idx]
MP = [[Bvecs[i][pw] for pw in pivots] for i in range(len(basis_idx))]
# invert exactly via Gauss-Jordan
n = len(MP)
A = [row[:] + [Fr(1) if i == j else Fr(0) for j in range(n)] for i, row in enumerate(MP)]
for col in range(n):
    piv = next(i for i in range(col, n) if A[i][col])
    A[[col, piv]] = A[[piv, col]]
    f = A[col][col]
    A[col] = [x/f for x in A[col]]
    for i in range(n):
        if i != col and A[i][col]:
            g = A[i][col]
            A[i] = [x - g*y for x, y in zip(A[i], A[col])]
MPinv = [row[n:] for row in A]
print("pivot matrix inverted; det check:", all(A[i][i] == 1 for i in range(n)))
pickle.dump({'basis_idx': basis_idx, 'pivots': pivots, 'MPinv': MPinv, 'Bn': Bn}, open('stage1.pkl','wb'))
