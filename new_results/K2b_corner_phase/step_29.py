
import sys; sys.path.insert(0,'research/08_theory/K2b_corner_phase/')
from k2utils import *
mp.dps = 80

def gam(a):  # gamma_+
    return mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4 - a/2))
def g1f(a): return -mp.mpf(3)/4 - 1j*mp.e**(1j*a)/(4*mp.sin(a))
def g2f(a):
    s = mp.sin(a); E = mp.e**(1j*a)
    return mp.mpf(29)/96 + E**2*(mp.mpf(1)/2 + 1/(32*s**2)) + (1j*E/s)*(E**2/4 - mp.mpf(3)/16)

def Jm(N, Z, m):
    # int_0^inf e^{-x} (1+Z/x)^{-N} (x+Z)^{m+3/2} dx/x
    xs = mp.sqrt(N*Z)
    f = lambda x: mp.e**(-x) * (1+Z/x)**(-N) * (x+Z)**(m+mp.mpf('1.5')) / x
    # integrate on log scale to capture saddle
    g = lambda t: f(mp.e**t)*mp.e**t
    return mp.quad(f, [0, xs], maxdegree=8) + mp.quad(f, [xs, mp.inf], maxdegree=8)

def Iplus(N, a, mmax=2):
    Z = 1j*a
    tot = 0
    gms = [mp.mpf(1), g1f(a), g2f(a)]
    for m in range(mmax+1):
        bf = mp.gamma(m+mp.mpf('1.5'))*mp.gamma(N-m-mp.mpf('1.5'))/mp.gamma(N)
        tot += gms[m]*bf*Jm(N,Z,m)
    return (-1j/mp.pi)*gam(a)*tot

def Wred(N, a):
    M = 2*mp.re(Iplus(N,a))
    A = mp.sqrt(2*a*mp.sin(a)); th = mp.sqrt(2*a*N)
    return -M*N*mp.e**th/A

a = mp.pi/2
for N in [100, 300]:
    print('Wred(N=%d) ='%N, mp.nstr(Wred(N,a),20))
    d = load('pi_2'); th = mp.sqrt(2*a*N)
    Wtab = -d[N-1]*N*mp.e**th/mp.sqrt(mp.pi)
    print('  table W   =', mp.nstr(mp.re(Wtab),20), ' (diff = image towers)')
