
import pickle, numpy as np
from fractions import Fraction as Fr
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve81.pkl','rb'))
li4v, prodv, Rv, bwords = D['li4v'], D['prodv'], D['Rv'], D['bwords']
S = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/stage1.pkl','rb'))
basis_idx, pivots, MPinv, Bn = S['basis_idx'], S['pivots'], S['MPinv'], S['Bn']
names = list(prodv.keys()); vecs = list(prodv.values())
Bvecs = [vecs[k] for k in basis_idx]
n = len(Bn)
nonpiv = [j for j in range(81) if j not in set(pivots)]

def Qimage(v):
    """v (81-dim Fr) -> coords on nonpivot words after subtracting product basis part"""
    vp = [v[pw] for pw in pivots]
    x = [sum(MPinv[i][j]*vp[j] for j in range(n)) for i in range(n)]
    resid = [v[j] - sum(x[i]*Bvecs[i][j] for i in range(n)) for j in range(81)]
    return [resid[j] for j in nonpiv], resid, x

# Stage 1: solve Q(R_E) = sum a_i Q(Li4_i) per exponent
li4n = list(li4v.keys())
Qli4 = []
for k in li4n:
    q, _, _ = Qimage(li4v[k])
    Qli4.append(q)
QR = {}
for E in (-1,0,1):
    q, _, _ = Qimage(Rv[E])
    QR[E] = q
# solve exactly: 18 eqs, 36 unknowns
import sympy as sp
a_sols = {}
for E in (-1,0,1):
    A = sp.Matrix(18, 36, lambda i,j: sp.Rational(Qli4[j][i]))
    b = sp.Matrix(18, 1, lambda i,j: sp.Rational(QR[E][i]))
    aug = A.row_join(b)
    rref_aug = aug.rref()[0]
    # identify pivot columns and free columns
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
    used = [(li4n[j], sol[j]) for j in range(36) if sol[j]]
    print(f"E={E}: {len(used)} Li4 terms used:", used)
    print("   free cols:", [li4n[j] for j in free])
pickle.dump({'a_sols': a_sols}, open('a_sols.pkl','wb'))
