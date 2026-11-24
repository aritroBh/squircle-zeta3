
import sys; sys.path.insert(0,'research/08_theory/K2b_corner_phase/')
from k2utils import *
mp.dps = 100

def fit_kappa(key, kmin=80, kmax=300, J=4):
    a = ALPHA[key]; d = load(key)
    A = mp.sqrt(2*a*mp.sin(a))
    ks = list(range(kmin,kmax+1))
    # build design: unknowns c=(Re k1, Im k1, ..., Re kJ, Im kJ)
    rows=[]; rhs=[]
    for k in ks:
        N = mp.mpf(k+1)
        th = mp.sqrt(2*a*N)
        W = -d[k]*N*mp.e**th/A
        y = W - mp.cos(th)
        row=[]
        for j in range(1,J+1):
            p = N**(-mp.mpf(j)/2)
            row += [p*mp.cos(th), -p*mp.sin(th)]
        rows.append(row); rhs.append(y)
    # least squares via normal equations at high precision
    n = 2*J
    G = mp.matrix(n,n); b = mp.matrix(n,1)
    for r in range(n):
        b[r] = sum(rows[i][r]*rhs[i] for i in range(len(ks)))
        for c in range(n):
            G[r,c] = sum(rows[i][r]*rows[i][c] for i in range(len(ks)))
    sol = mp.lu_solve(G,b)
    kappas = [mp.mpc(sol[2*j-2], sol[2*j-1]) for j in range(1,J+1)]
    return kappas

k1_true = ((24+12*mp.pi+mp.pi**2) + 1j*(24-12*mp.pi+mp.pi**2))/(48*mp.sqrt(mp.pi))
print('kappa1(pi/2) known:', mp.nstr(k1_true, 20))
for J in [3,4,5]:
    ks_ = fit_kappa('pi_2', 80, 300, J)
    print('pi_2 J=%d:'%J, ' '.join(mp.nstr(k,15) for k in ks_[:2]))
