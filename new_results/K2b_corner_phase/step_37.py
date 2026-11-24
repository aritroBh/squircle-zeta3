
import sys; sys.path.insert(0,'research/08_theory/K2b_corner_phase/')
from k2utils import *
mp.dps = 100
def Jm(N, Z, m):
    xs = mp.sqrt(N*Z)
    f = lambda x: mp.e**(-x) * (1+Z/x)**(-N) * (x+Z)**(m+mp.mpf('1.5')) / x
    return mp.quad(f, [0, xs], maxdegree=9) + mp.quad(f, [xs, mp.inf], maxdegree=9)
def bracket_true(N, Z, m):
    xs = (-Z + mp.sqrt(Z**2 + 4*N*Z))/2
    phis = xs + N*mp.log(1+Z/xs)
    phi2 = N*(1/xs**2 - 1/(xs+Z)**2)
    h = (xs+Z)**(m+mp.mpf('1.5'))/xs
    G = mp.sqrt(2*mp.pi/phi2)
    return Jm(N,Z,m) / (mp.e**(-phis)*h*G)
Z = 1j*mp.pi/2
for m in [0,2]:
    for N in [2000, 4000, 8000]:
        print('m=%d N=%5d bracket=%s'%(m,N,mp.nstr(bracket_true(mp.mpf(N),Z,m),20)))
# expected b1 m=0: 1/(2 sqZ) = 0.282-0.282i ; sig=1/sqrt(4000)=0.0158 -> b1 sig ~ 0.00446-0.00446i
