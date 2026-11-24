
import sys; sys.path.insert(0,'research/08_theory/K2b_corner_phase/')
from k2utils import *
mp.dps = 80
def Jm(N, Z, m):
    xs = mp.sqrt(N*Z)
    f = lambda x: mp.e**(-x) * (1+Z/x)**(-N) * (x+Z)**(m+mp.mpf('1.5')) / x
    return mp.quad(f, [0, xs], maxdegree=9) + mp.quad(f, [xs, mp.inf], maxdegree=9)
def bracket_true(N, Z, m):
    # exact saddle quantities
    xs = (-Z + mp.sqrt(Z**2 + 4*N*Z))/2
    phis = xs + N*mp.log(1+Z/xs)
    phi2 = N*(1/xs**2 - 1/(xs+Z)**2)
    h = (xs+Z)**(m+mp.mpf('1.5'))/xs
    G = mp.sqrt(2*mp.pi/phi2)
    return Jm(N,Z,m) / (mp.e**(-phis)*h*G)
Z = 1j*mp.pi/2
for m in [0,1,2]:
    Ns = [3000, 5000, 8000, 12000]
    bs = [bracket_true(mp.mpf(N),Z,m) for N in Ns]
    # fit 1 + b1 sig + b2 sig^2 + b3 sig^3
    rows=[[ (1/mp.sqrt(N))**j for j in range(1,4)] for N in Ns]
    rhs=[bs[i]-1 for i in range(len(Ns))]
    G = mp.matrix(3,3); bb = mp.matrix(3,1)
    for r in range(3):
        bb[r] = sum(mp.conj(rows[i][r])*rhs[i] for i in range(len(Ns)))
        for c in range(3):
            G[r,c] = sum(mp.conj(rows[i][r])*rows[i][c] for i in range(len(Ns)))
    sol = mp.lu_solve(G,bb)
    print('m=%d true bracket: b1=%s b2=%s b3=%s'%(m, mp.nstr(sol[0],16), mp.nstr(sol[1],16), mp.nstr(sol[2],12)))
sqZ = mp.sqrt(Z)
print('my b1 m=0: 1/(2 sqZ) =', mp.nstr(1/(2*sqZ),16))
print('my b2 m=0: 2835/(512 Z) - 5/8 =', mp.nstr(mp.mpf(2835)/(512*Z)-mp.mpf(5)/8,16))
print('my b2 m=1: 699/(512Z) - 7/4 =', mp.nstr(mp.mpf(699)/(512*Z)-mp.mpf(7)/4,16))
print('my b2 m=2: -285/(512Z) -27/8 =', mp.nstr(-mp.mpf(285)/(512*Z)-mp.mpf(27)/8,16))
