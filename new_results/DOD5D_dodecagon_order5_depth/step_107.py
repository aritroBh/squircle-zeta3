import sys, pickle, time
sys.path.insert(0, ".")
from fractions import Fraction as F
import flint
import obstruction32 as OB
from qs3_symbol import Q3V

D = pickle.load(open("dod5d_symbols.pkl","rb"))
vecdata = pickle.load(open("mpl_cobrackets.pkl","rb"))
def height(v):
    h = 0
    for vs in v:
        vs = Q3V(vs) if not isinstance(vs, tuple) else vs
        for x in (vs[0], vs[1]):
            h = max(h, abs(x.numerator), x.denominator)
    return h
vd = sorted(vecdata, key=lambda t: height(t[3]))
print("height range:", height(vd[0][3]), "to", height(vd[-1][3]))
# keep the 1200 lowest-height vectors
vd = vd[:1200]
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
t0 = time.time()
R, rank = A.rref()
print("rref rank", rank, "%.0fs" % (time.time()-t0))
inconsistent = any(R[r, ncols] != 0 and all(R[r, j] == 0 for j in range(ncols)) for r in range(48))
print("inconsistent:", inconsistent)
pivlist = []
for r in range(rank):
    for j in range(ncols):
        if R[r, j] != 0: pivlist.append(j); break
x = [flint.fmpq(0)]*ncols
for r, pc in enumerate(pivlist):
    if pc < ncols: x[pc] = R[r, ncols]
nz = [(j, x[j]) for j in range(ncols) if x[j] != 0]
print("nonzero coefficients:", len(nz))
res = [sum(A[i, j]*x[j] for j in range(ncols)) - A[i, ncols] for i in range(48)]
print("residual zero:", all(r == 0 for r in res))
terms = {}
for j, xj in nz:
    terms.setdefault(j//2, [F(0), F(0)])[j % 2] = F(xj)
def fmtq3(p, q):
    if q == 0: return str(p)
    if p == 0: return f"{q}r3"
    return f"{p}{q:+}r3"
for i, (p, q) in sorted(terms.items()):
    (ab, xx, yy, _) = vd[i]
    print(f"  ({fmtq3(p,q)}) * Li_{ab[0]},{ab[1]}( ({fmtq3(F(xx[0]),F(xx[1]))}), ({fmtq3(F(yy[0]),F(yy[1]))}) )")
pickle.dump([(i, str(terms[i][0]), str(terms[i][1])) for i in sorted(terms)], open("li32_solution_clean.pkl","wb"))
# save the vector subset used
pickle.dump(vd, open("mpl_cobrackets_low.pkl","wb"))
