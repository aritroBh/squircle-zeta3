
import sys; sys.path.insert(0,'research/08_theory/K2b_corner_phase/')
from k2utils import *
mp.dps = 100
def Jm(N, Z, m):
    xs = mp.sqrt(N*Z)
    dl = 12*N**mp.mpf('0.25')
    f = lambda x: mp.e**(-x) * (1+Z/x)**(-N) * (x+Z)**(m+mp.mpf('1.5')) / x
    return (mp.quad(f, [0, xs-dl], maxdegree=9) + mp.quad(f, [xs-dl, xs, xs+dl], maxdegree=9)
            + mp.quad(f, [xs+dl, mp.inf], maxdegree=9))
def bracket_true(N, Z, m):
    xs = (-Z + mp.sqrt(Z**2 + 4*N*Z))/2
    phis = xs + N*mp.log(1+Z/xs)
    phi2 = N*(1/xs**2 - 1/(xs+Z)**2)
    h = (xs+Z)**(m+mp.mpf('1.5'))/xs
    G = mp.sqrt(2*mp.pi/phi2)
    return Jm(N,Z,m) / (mp.e**(-phis)*h*G)
Z = 1j*mp.pi/2
for m in [0,1,2]:
    Ns = [2000, 4000, 8000, 16000]
    bs = [bracket_true(mp.mpf(N),Z,m) for N in Ns]
    rows=[[ (1/mp.sqrt(N))**j for j in range(1,5)] for N in Ns]
    rhs=[bs[i]-1 for i in range(len(Ns))]
    G = mp.matrix(4,4); bb = mp.matrix(4,1)
    for r in range(4):
        bb[r] = sum(mp.conj(rows[i][r])*rhs[i] for i in range(len(Ns)))
        for c in range(4):
            G[r,c] = sum(mp.conj(rows[i][r])*rows[i][c] for i in range(len(Ns)))
    sol = mp.lu_solve(G,bb)
    print('m=%d  b1=%s'%(m, mp.nstr(sol[0],20)))
    print('      b2=%s'%mp.nstr(sol[1],20))
    print('      b3=%s  b4=%s'%(mp.nstr(sol[2],14), mp.nstr(sol[3],12)))
