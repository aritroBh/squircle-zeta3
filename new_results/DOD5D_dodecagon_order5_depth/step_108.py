import sys, pickle
sys.path.insert(0, ".")
from fractions import Fraction as F
import flint
import obstruction32 as OB
from qs3_symbol import Q3V

D = pickle.load(open("dod5d_symbols.pkl","rb"))
vd = pickle.load(open("mpl_cobrackets_low.pkl","rb"))
L3 = OB.lie3_basis(3); L2 = OB.lie2_basis(3)
Psi = OB.psi_tensor(D['total'])
obs = [Q3V(OB.pair(Psi, A, B)) for (_, A) in L3 for (_, B) in L2]
N = len(vd); ncols = 2*N
rows = []
for s in range(24):
    rowp, rowq = [], []
    for i, (_, _, _, v) in enumerate(vd):
        vs = Q3V(v[s]) if not isinstance(v[s], tuple) else v[s]
        rowp += [F(vs[0]), 3*F(vs[1])]
        rowq += [F(vs[1]), F(vs[0])]
    rows.append(rowp + [F(obs[s][0])]); rows.append(rowq + [F(obs[s][1])])
A = flint.fmpq_mat([[flint.fmpq(str(x)) for x in row] for row in rows])
R, rank = A.rref()
pivlist = []
for r in range(rank):
    for j in range(ncols):
        if R[r, j] != 0: pivlist.append(j); break
x = [flint.fmpq(0)]*ncols
for r, pc in enumerate(pivlist):
    if pc < ncols: x[pc] = R[r, ncols]
terms = {}
for j in range(ncols):
    if x[j] != 0:
        terms.setdefault(j//2, [F(0), F(0)])[j % 2] = F(str(x[j]))
def fmtq3(p, q):
    if q == 0: return str(p)
    if p == 0: return f"{q}r3"
    return f"{p}{q:+}r3"
lines = []
for i, (p, q) in sorted(terms.items()):
    (ab, xx, yy, _) = vd[i]
    lines.append(f"  ({fmtq3(p,q)}) * Li_{ab[0]},{ab[1]}( ({fmtq3(F(xx[0]),F(xx[1]))}), ({fmtq3(F(yy[0]),F(yy[1]))}) )")
print("\n".join(lines))
print("terms:", len(terms))
# verify again explicitly
res = [sum(A[i, j]*x[j] for j in range(ncols)) - A[i, ncols] for i in range(48)]
print("residual zero:", all(r == 0 for r in res))
pickle.dump(lines, open("li32_solution_lines.pkl","wb"))
