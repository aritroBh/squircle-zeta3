
import sys; sys.path.insert(0,'research/08_theory/K2b_corner_phase/')
from k2utils import *
mp.dps = 80
def gam(a): return mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4 - a/2))
def g1f(a): return -mp.mpf(3)/4 - 1j*mp.e**(1j*a)/(4*mp.sin(a))
def g2f(a):
    s = mp.sin(a); E = mp.e**(1j*a)
    return mp.mpf(29)/96 + E**2*(mp.mpf(1)/2 + 1/(32*s**2)) + (1j*E/s)*(E**2/4 - mp.mpf(3)/16)
def Jm(N, Z, m):
    xs = mp.sqrt(N*Z)
    f = lambda x: mp.e**(-x) * (1+Z/x)**(-N) * (x+Z)**(m+mp.mpf('1.5')) / x
    return mp.quad(f, [0, xs], maxdegree=9) + mp.quad(f, [xs, mp.inf], maxdegree=9)
def Iplus(N, a, mmax=2):
    Z = 1j*a; tot = 0
    gms = [mp.mpf(1), g1f(a), g2f(a)]
    for m in range(mmax+1):
        bf = mp.gamma(m+mp.mpf('1.5'))*mp.gamma(N-m-mp.mpf('1.5'))/mp.gamma(N)
        tot += gms[m]*bf*Jm(N,Z,m)
    return (+1j/mp.pi)*gam(a)*tot
def Pphasor(N, a):
    A = mp.sqrt(2*a*mp.sin(a)); th = mp.sqrt(2*a*N)
    return -(2/A)*N*mp.e**(th*(1+1j))*Iplus(N,a)

a = mp.pi/2
Ns = [2000,3000,4000,5000,7000,10000,15000,20000,30000,50000]
Ps = []
for N in Ns:
    Ps.append(Pphasor(N,a))
    print('N=%6d P=%s'%(N, mp.nstr(Ps[-1],18)))
# least squares P = sum_{j=0..5} kap_j sig^j
J = 5
rows=[]
for N,P in zip(Ns,Ps):
    sg = 1/mp.sqrt(N)
    rows.append([sg**j for j in range(J+1)])
# complex LS: solve min ||rows kap - P||
G = mp.matrix(J+1,J+1); b = mp.matrix(J+1,1)
for r in range(J+1):
    b[r] = sum(mp.conj(rows[i][r])*Ps[i] for i in range(len(Ns)))
    for c in range(J+1):
        G[r,c] = sum(mp.conj(rows[i][r])*rows[i][c] for i in range(len(Ns)))
sol = mp.lu_solve(G,b)
for j in range(J+1):
    print('kappa_%d = %s'%(j, mp.nstr(sol[j],20)))
print('my kappa1 =', mp.nstr(mp.mpc(0.84121508776758439708,-0.04501183768517361657),20))
