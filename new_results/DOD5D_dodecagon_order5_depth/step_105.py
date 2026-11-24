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
A = flint.fmpq_mat([[flint.fmpq(str(x)) for x in row] for row in rows])
t0 = time.time()
R, rank = A.rref()
print("rref: rank", rank, "%.0fs" % (time.time()-t0))
inconsistent = False
for r in range(nrows):
    if R[r, ncols] != 0 and all(R[r, j] == 0 for j in range(ncols)):
        inconsistent = True
print("system inconsistent:", inconsistent)
pivlist = []
for r in range(rank):
    for j in range(ncols):
        if R[r, j] != 0:
            pivlist.append(j); break
x = [flint.fmpq(0)]*ncols
for r, pc in enumerate(pivlist):
    if pc < ncols:
        x[pc] = R[r, ncols]
nz = [(j, x[j]) for j in range(ncols) if x[j] != 0]
print("nonzero coefficients:", len(nz))
t0 = time.time()
res = [sum(A[i, j]*x[j] for j in range(ncols)) - A[i, ncols] for i in range(nrows)]
print("residual zero:", all(r == 0 for r in res), "%.0fs" % (time.time()-t0))
pickle.dump([(j, str(xj)) for j, xj in nz], open("li32_solution_raw.pkl","wb"))
