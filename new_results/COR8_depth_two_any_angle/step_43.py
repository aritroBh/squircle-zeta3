
import pickle
from fractions import Fraction as Fr
import sympy as sp
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve81.pkl','rb'))
li4v, prodv, Rv = D['li4v'], D['prodv'], D['Rv']
S = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/stage1.pkl','rb'))
basis_idx, pivots, MPinv = S['basis_idx'], S['pivots'], S['MPinv']
names = list(prodv.keys()); vecs = list(prodv.values())
Bvecs = [vecs[k] for k in basis_idx]
n = len(Bn) if False else len(basis_idx)
nonpiv = [j for j in range(81) if j not in set(pivots)]

def Qimage(v):
    vp = [v[pw] for pw in pivots]
    x = [sum(MPinv[i][j]*vp[j] for j in range(n)) for i in range(n)]
    resid = [v[j] - sum(x[i]*Bvecs[i][j] for i in range(n)) for j in range(81)]
    return [resid[j] for j in nonpiv], resid, x

li4n = list(li4v.keys())
Qli4 = [Qimage(li4v[k])[0] for k in li4n]
a_sols = {}
for E in (-1,0,1):
    qR, _, _ = Qimage(Rv[E])
    A = sp.Matrix(18, 36, lambda i,j: sp.Rational(Qli4[j][i]))
    b = sp.Matrix(18, 1, lambda i,j: sp.Rational(qR[i]))
    rref_aug = A.row_join(b).rref()[0]
    piv_cols = []
    for i in range(18):
        for j in range(36):
            if rref_aug[i,j] != 0:
                piv_cols.append(j); break
    free = [j for j in range(36) if j not in piv_cols]
    sol = [Fr(0)]*36
    for i, j in enumerate(piv_cols):
        sol[j] = Fr(rref_aug[i, 36])
    a_sols[E] = sol
    print(f"E={E}: used:", [(li4n[j], str(sol[j])) for j in range(36) if sol[j]])
    print("   free:", [li4n[j] for j in free])
pickle.dump({'a_sols': a_sols, 'li4n': li4n}, open('a_sols.pkl','wb'))
