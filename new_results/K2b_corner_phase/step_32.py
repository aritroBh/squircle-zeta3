
import sys; sys.path.insert(0,'research/08_theory/K2b_corner_phase/')
from k2utils import *
mp.dps = 80
def gam(a): return mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4 - a/2))
def g1f(a): return -mp.mpf(3)/4 - 1j*mp.e**(1j*a)/(4*mp.sin(a))
def g2f(a):
    s = mp.sin(a); E = mp.e**(1j*a)
    return mp.mpf(29)/96 + E**2*(mp.mpf(1)/2 + 1/(32*s**2)) + (1j*E/s)*(E**2/4 - mp.mpf(3)/16)
def Jm(N, Z, m):
    sqNZ = mp.sqrt(N*Z); y0 = 1 - mp.sqrt(Z)/(2*mp.sqrt(N))  # approx saddle location
    f = lambda x: mp.e**(-x) * (1+Z/x)**(-N) * (x+Z)**(m+mp.mpf('1.5')) / x
    xs = sqNZ*y0
    w = 8*(N)**mp.mpf('0.25')
    return mp.quad(f, [0, xs], maxdegree=8) + mp.quad(f, [xs-w, xs, xs+w], maxdegree=8) + mp.quad(f, [xs, mp.inf], maxdegree=8)
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
for N in [1000, 10000]:
    P = Pphasor(N,a)
    print('N=%6d  P=%s' % (N, mp.nstr(P, 25)))
    print('   expected ~ 1 + k1/sig: k1/sig =', mp.nstr(mp.mpc(0.84121508776758439708,-0.04501183768517361657)/mp.sqrt(N),18))
