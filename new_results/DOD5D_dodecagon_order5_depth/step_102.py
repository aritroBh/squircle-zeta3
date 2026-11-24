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
print("obs entries:", len(obs), "nonzero:", sum(1 for v in obs if v != 0))

N = len(vecdata)
nrows, ncols = 48, 2*N
t0 = time.time()
rows = []
for s in range(24):
    rp = flint.fmpq(0); 
    rowp = []
    rowq = []
    for i, (_, _, _, v) in enumerate(vecdata):
        vp, vq = F(v[s][0]), F(v[s][1])
        rowp += [vp, 3*vq]
        rowq += [vq, vp]
    rows.append(rowp); rows.append(rowq)
b = []
for s in range(24):
    b += [F(obs[s][0]), F(obs[s][1])]
print("built rows %.0fs" % (time.time()-t0))
A = flint.fmpq_mat([[flint.fmpq(str(x)) for x in row] for row in rows])
bb = flint.fmpq_mat([[flint.fmpq(str(x))] for x in b])
aug = A.concat(bb)
print("aug", aug.nrows(), aug.ncols(), "%.0fs" % (time.time()-t0))
t0 = time.time()
R, rank, pivots = aug.rref()
print("rref: rank", rank, "%.0fs" % (time.time()-t0))
# consistency: any row with zero A-part but nonzero last col?
inconsistent = False
for r in range(rank, nrows):
    if R[r, ncols] != 0:
        inconsistent = True
print("system inconsistent (obstruction NOT in span):", inconsistent)
pivlist = sorted(pivots)
print("pivot cols:", pivlist[:10], "...")
# particular solution: free vars 0
sol = {}
for r, pc in enumerate(pivlist):
    if pc < ncols:
        sol[pc] = R[r, ncols]
x = [sol.get(j, flint.fmpq(0)) for j in range(ncols)]
# verify A x = b
res = [sum(A[i, j]*x[j] for j in range(ncols)) - bb[i, 0] for i in range(nrows)]
print("residual zero:", all(r == 0 for r in res))
nz = [(j, x[j]) for j in range(ncols) if x[j] != 0]
print("nonzero coefficients:", len(nz))
pickle.dump([(j, str(xj)) for j, xj in nz], open("li32_solution_raw.pkl","wb"))
