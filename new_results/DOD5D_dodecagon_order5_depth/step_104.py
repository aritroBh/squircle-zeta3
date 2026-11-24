import sys, pickle, time
sys.path.insert(0, ".")
from fractions import Fraction as F
import flint
import obstruction32 as OB
from qs3_symbol import Q3V

D = pickle.load(open("dod5d_symbols.pkl","rb"))
vecdata = pickle.load(open("mpl_cobrackets.pkl","rb"))
L3 = OB.lie3_basis(3); L2 = OB.lie2_basis(3)
Psi = OB.psi_tensor(D['total'])
obs = [Q3V(OB.pair(Psi, A, B)) for (_, A) in L3 for (_, B) in L2]
N = len(vecdata)
nrows, ncols = 48, 2*N
rows = []
for s in range(24):
    rowp, rowq = [], []
    for i, (_, _, _, v) in enumerate(vecdata):
        vs = Q3V(v[s]) if not isinstance(v[s], tuple) else v[s]
        vp, vq = F(vs[0]), F(vs[1])
        rowp += [vp, 3*vq]
        rowq += [vq, vp]
    rows.append(rowp + [F(obs[s][0])]); rows.append(rowq + [F(obs[s][1])])
t0 = time.time()
A = flint.fmpq_mat([[flint.fmpq(str(x)) for x in row] for row in rows])
print("built aug", A.nrows(), A.ncols(), "%.0fs" % (time.time()-t0))
t0 = time.time()
R, rank, pivots = A.rref()
print("rref: rank", rank, "%.0fs" % (time.time()-t0))
inconsistent = any(R[r, ncols] != 0 and all(R[r, j] == 0 for j in range(ncols)) for r in range(nrows))
print("system inconsistent (obstruction NOT in span):", inconsistent)
pivlist = sorted(pivots)
sol = {}
for r, pc in enumerate(pivlist):
    if pc < ncols:
        sol[pc] = R[r, ncols]
x = [sol.get(j, flint.fmpq(0)) for j in range(ncols)]
nz = [(j, x[j]) for j in range(ncols) if x[j] != 0]
print("nonzero coefficients:", len(nz))
pickle.dump([(j, str(xj)) for j, xj in nz], open("li32_solution_raw.pkl","wb"))
# verify residual
t0 = time.time()
res = [sum(A[i, j]*x[j] for j in range(ncols)) - A[i, ncols] for i in range(nrows)]
print("residual zero:", all(r == 0 for r in res), "%.0fs" % (time.time()-t0))
