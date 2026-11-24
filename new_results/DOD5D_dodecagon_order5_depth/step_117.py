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
terms = []
for j in range(ncols):
    if x[j] != 0:
        i = j//2
        ab, xx, yy, _ = vd[i]
        terms.append((ab, xx, yy, j % 2, F(str(x[j]))))
# consolidate
tmap = {}
for ab, xx, yy, pq, c in terms:
    key = (ab, xx, yy)
    tmap.setdefault(key, [F(0), F(0)])[pq] = c
sol_terms = [(ab, xx, yy, Q3V(c[0], c[1])) for (ab, xx, yy), c in tmap.items()]
pickle.dump(sol_terms, open("li32_solution.pkl","wb"))
print("saved li32_solution.pkl:", len(sol_terms), "terms")
# independent re-verification of the certificate from the saved file
combo = [Q3V(0)]*24
for ab, xx, yy, c in sol_terms:
    for (ab2, xx2, yy2, v) in vd:
        if ab2 == ab and xx2 == xx and yy2 == yy:
            for s in range(24):
                combo[s] = combo[s] + c*Q3V(v[s])
            break
resid = [combo[s] - obs[s] for s in range(24)]
print("certificate recomputed from file: residual zero:", all(r == 0 for r in resid))
