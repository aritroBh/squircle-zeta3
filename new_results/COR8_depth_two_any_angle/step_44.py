
import pickle
from fractions import Fraction as Fr
import sympy as sp
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve81.pkl','rb'))
li4v, prodv, Rv = D['li4v'], D['prodv'], D['Rv']
S = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/stage1.pkl','rb'))
basis_idx, pivots, MPinv = S['basis_idx'], S['pivots'], S['MPinv']
names = list(prodv.keys()); vecs = list(prodv.values())
Bvecs = [vecs[k] for k in basis_idx]
n = len(basis_idx)
nonpiv = [j for j in range(81) if j not in set(pivots)]
def Qimage(v):
    vp = [v[pw] for pw in pivots]
    x = [sum(MPinv[i][j]*vp[j] for j in range(n)) for i in range(n)]
    resid = [v[j] - sum(x[i]*Bvecs[i][j] for i in range(n)) for j in range(81)]
    return [resid[j] for j in nonpiv]
li4n = list(li4v.keys())
Qli4 = {k: Qimage(li4v[k]) for k in li4n}
nice = ['c', '-c', '1 - c', '1/(c + 1)', 'c + 1', '(1 - c)/(c + 1)', '(c - 1)/(c + 1)', '(1 - c)*(c + 1)',
        '-c/(1 - c)', 'c/(c + 1)', 'c**2', '-c**2/((1 - c)*(c + 1))', '2/(1 - c)', '2/(c + 1)',
        '-2*c/(1 - c)', '2*c/(c + 1)', '-4*c/(1 - c)**2', '4*c/(c + 1)**2', '-1/c', '-(c + 1)/(1 - c)',
        '1/(1 - c)', '1/c', '(c + 1)/c', '-(1 - c)/c', '(1 - c)**2/(c + 1)**2']
nice = [s for s in nice if s in Qli4]
print(len(nice), "nice args")
a_sols = {}
for E in (-1,0,1):
    A = sp.Matrix(18, len(nice), lambda i,j: sp.Rational(Qli4[nice[j]][i]))
    qR = Qimage(Rv[E])
    b = sp.Matrix(18, 1, lambda i,j: sp.Rational(qR[i]))
    M = A.row_join(b)
    r = M.rref()[0]
    # check consistency and read particular solution
    inconsistent = any(all(r[i,j]==0 for j in range(len(nice))) and r[i,len(nice)]!=0 for i in range(18))
    print(f"E={E}: inconsistent={inconsistent}, rank A = {A.rank()}")
    if inconsistent:
        continue
    piv_cols = []
    for i in range(18):
        for j in range(len(nice)):
            if r[i,j] != 0:
                piv_cols.append(j); break
    sol = [Fr(0)]*len(nice)
    for i, j in enumerate(piv_cols):
        sol[j] = Fr(r[i, len(nice)])
    a_sols[E] = sol
    print("   used:", [(nice[j], str(sol[j])) for j in range(len(nice)) if sol[j]])
pickle.dump({'a_sols': a_sols, 'nice': nice}, open('a_sols.pkl','wb'))
