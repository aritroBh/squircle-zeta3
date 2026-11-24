
import sys; sys.path.insert(0,'research/08_theory/K2b_corner_phase/')
from k2utils import *
mp.dps = 80
def g1f(a): return -mp.mpf(3)/4 - 1j*mp.e**(1j*a)/(4*mp.sin(a))
def g2f(a):
    s = mp.sin(a); E = mp.e**(1j*a)
    return mp.mpf(29)/96 + E**2*(mp.mpf(1)/2 + 1/(32*s**2)) + (1j*E/s)*(E**2/4 - mp.mpf(3)/16)
def beta1(a):
    Z = 1j*a; sqZ = mp.sqrt(Z)
    return 1/(2*sqZ) + sqZ/8 - Z**mp.mpf('1.5')/12 - 3j*mp.e**(1j*a)*sqZ/(8*mp.sin(a))
def beta2(a):
    Z = 1j*a; g1 = g1f(a); g2 = g2f(a)
    return (Z**3/288 - Z**2*g1/8 - 5*Z**2/48 + 21*Z*g1/8 + 15*Z*g2/4 + 65*Z/96 + 9*g1/4 + mp.mpf(15)/8 + mp.mpf(2835)/(512*Z))
def kappa1(a): return mp.conj(beta1(a))
def kappa2(a): return mp.conj(beta2(a))

print('%-7s %5s | %-14s %-14s | %-12s %-12s | %-12s'%('alpha','k','rho0 (lead)','rho0*N^.5','rho1 (2term)','rho1*N','rho2 (3term)*N^1.5'))
for key in ['pi_3','pi_2','2pi_3','3pi_4','5pi_6']:
    a = ALPHA[key]; d = load(key); A = mp.sqrt(2*a*mp.sin(a))
    k1 = kappa1(a); k2 = kappa2(a)
    for k in [100,150,200,250,300]:
        N = mp.mpf(k+1); th = mp.sqrt(2*a*N); sqN = mp.sqrt(N)
        W = -d[k]*N*mp.e**th/A
        c, s = mp.cos(th), mp.sin(th)
        def ReEi(x): return mp.re(x)*c - mp.im(x)*s
        pred0 = c
        pred1 = c + ReEi(k1)/sqN
        pred2 = pred1 + ReEi(k2)/N
        rho0 = W - pred0
        rho1 = W - pred1
        rho2 = W - pred2
        print('%-7s %5d | %-14s %-14s | %-12s %-12s | %-12s'%(key, k,
              mp.nstr(mp.re(rho0),8), mp.nstr(mp.re(rho0*sqN),8),
              mp.nstr(mp.re(rho1),8), mp.nstr(mp.re(rho1*N),8), mp.nstr(mp.re(rho2*N**mp.mpf('1.5')),8)))
