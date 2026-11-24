
import sys; sys.path.insert(0,'research/08_theory/K2b_corner_phase/')
from k2utils import *
mp.dps = 80
Ns = [2000,3000,4000,5000,7000,10000]
Ps = [(1.01923377454113546 + 0.00113376997429773435j),
(1.01564124186477023 + 0.000904064451203547315j),
(1.01351306737739058 + 0.000772276477237001862j),
(1.01206649797325144 + 0.000684424156525533428j),
(1.01017589509898513 + 0.00057162741706328805j),
(1.00849721511378789 + 0.000473322815545924496j)]
Ps = [mp.mpc(p) for p in Ps]
def fit(J):
    rows=[[ (1/mp.sqrt(N))**j for j in range(1,J+1)] for N in Ns]
    rhs=[Ps[i]-1 for i in range(len(Ns))]
    G = mp.matrix(J,J); b = mp.matrix(J,1)
    for r in range(J):
        b[r] = sum(mp.conj(rows[i][r])*rhs[i] for i in range(len(Ns)))
        for c in range(J):
            G[r,c] = sum(mp.conj(rows[i][r])*rows[i][c] for i in range(len(Ns)))
    return mp.lu_solve(G,b)
for J in [2,3,4]:
    sol = fit(J)
    print('J=%d: k1=%s k2=%s'%(J, mp.nstr(sol[0],18), mp.nstr(sol[1],18)))
print('target k1 = 0.84121508776758439708 - 0.04501183768517361657j')
print('my (buggy) k2 = 0.85280837917801415228 + 3.3094164309273583075j')
