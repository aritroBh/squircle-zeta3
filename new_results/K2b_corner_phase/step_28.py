
import sys; sys.path.insert(0,'research/08_theory/K2b_corner_phase/')
from k2utils import *
mp.dps = 80
def g1f(a): return -mp.mpf(3)/4 - 1j*mp.e**(1j*a)/(4*mp.sin(a))
def g2f(a):
    s = mp.sin(a); E = mp.e**(1j*a)
    return mp.mpf(29)/96 + E**2*(mp.mpf(1)/2 + 1/(32*s**2)) + (1j*E/s)*(E**2/4 - mp.mpf(3)/16)
def beta2(a):
    Z = 1j*a; g1 = g1f(a); g2 = g2f(a)
    return (Z**3/288 - Z**2*g1/8 - 5*Z**2/48 + 21*Z*g1/8 + 15*Z*g2/4 + 65*Z/96 + 9*g1/4 + mp.mpf(15)/8 + mp.mpf(2835)/(512*Z))
def kappa2(a): return mp.conj(beta2(a))
a = mp.pi/2; d = load('pi_2'); A = mp.sqrt(mp.pi)
cr = 1/(2*mp.sqrt(mp.pi)) + mp.sqrt(mp.pi)/4 + mp.pi**mp.mpf('1.5')/48
ci = -1/(2*mp.sqrt(mp.pi)) + mp.sqrt(mp.pi)/4 - mp.pi**mp.mpf('1.5')/48
k1 = mp.mpc(cr, -ci)
k2 = kappa2(a)
print('my kappa2(pi/2) =', mp.nstr(k2,20))
print('k    observed rho1*N      predicted Re[e^{i th} k2]    difference')
for k in range(95, 301, 15):
    N = mp.mpf(k+1); th = mp.sqrt(mp.pi*N); sqN = mp.sqrt(N)
    W = -d[k]*N*mp.e**th/A
    c, s = mp.cos(th), mp.sin(th)
    rho1 = W - (c + (mp.re(k1)*c - mp.im(k1)*s)/sqN)
    obs = rho1*N
    pred = mp.re(k2)*c - mp.im(k2)*s
    print('%3d  %-18s %-18s %-18s' % (k, mp.nstr(mp.re(obs),12), mp.nstr(pred,12), mp.nstr(mp.re(obs-pred),12)))
