
import pickle, numpy as np
from fractions import Fraction as Fr
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve81.pkl','rb'))
li4v, prodv, Rv, bwords = D['li4v'], D['prodv'], D['Rv'], D['bwords']
names = list(prodv.keys()); vecs = list(prodv.values())
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
basis_idx, B_rows, pivots = [], [], []
for k in order:
    v = np.array([x % p for x in vecs[k]], dtype=np.int64)
    for br, pw in zip(B_rows, pivots):
        if v[pw]: v = (v - int(v[pw])*br) % p
    nz = np.nonzero(v)[0]
    if len(nz):
        pw = int(nz[0]); B_rows.append((v*pow(int(v[pw]), p-2, p)) % p); pivots.append(pw); basis_idx.append(k)
n = len(basis_idx)
print("product basis:", n)
Bvecs = [vecs[k] for k in basis_idx]
MP = [[Bvecs[i][pw] for pw in pivots] for i in range(n)]
A = [row[:] + [Fr(1) if i == j else Fr(0) for j in range(n)] for i, row in enumerate(MP)]
for col in range(n):
    piv = next(i for i in range(col, n) if A[i][col])
    A[col], A[piv] = A[piv], A[col]
    f = A[col][col]; A[col] = [x/f for x in A[col]]
    for i in range(n):
        if i != col and A[i][col]:
            g = A[i][col]; A[i] = [x - g*y for x, y in zip(A[i], A[col])]
MPinv = [row[n:] for row in A]
pickle.dump({'basis_idx': basis_idx, 'pivots': pivots, 'MPinv': MPinv}, open('stage1.pkl','wb'))
print("saved stage1.pkl")
