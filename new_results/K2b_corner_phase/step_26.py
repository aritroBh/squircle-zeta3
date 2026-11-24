
import sys; sys.path.insert(0,'research/08_theory/K2b_corner_phase/')
from k2utils import *
mp.dps = 80

def g1f(a):
    return -mp.mpf(3)/4 - 1j*mp.e**(1j*a)/(4*mp.sin(a))
def g2f(a):
    s = mp.sin(a); E = mp.e**(1j*a)
    return mp.mpf(29)/96 + E**2*(mp.mpf(1)/2 + 1/(32*s**2)) + (1j*E/s)*(E**2/4 - mp.mpf(3)/16)

def beta1(a):
    Z = 1j*a; sqZ = mp.sqrt(Z)
    return 1/(2*sqZ) + sqZ/8 - Z**mp.mpf('1.5')/12 - 3j*mp.e**(1j*a)*sqZ/(8*mp.sin(a))
def beta2(a):
    Z = 1j*a; g1 = g1f(a); g2 = g2f(a)
    return (Z**3/288 - Z**2*g1/8 - 5*Z**2/48 + 21*Z*g1/8 + 15*Z*g2/4
            + 65*Z/96 + 9*g1/4 + mp.mpf(15)/8 + mp.mpf(2835)/(512*Z))

def kappa1(a): return mp.conj(beta1(a))
def kappa2(a): return mp.conj(beta2(a))

# calibration at pi/2
k1_true = ((24+12*mp.pi+mp.pi**2) + 1j*(24-12*mp.pi+mp.pi**2))/(48*mp.sqrt(mp.pi))
print('kappa1(pi/2) mine :', mp.nstr(kappa1(mp.pi/2), 25))
print('kappa1(pi/2) D13b :', mp.nstr(k1_true, 25))
print('diff:', mp.nstr(abs(kappa1(mp.pi/2)-k1_true),5))
print('kappa2(pi/2) mine :', mp.nstr(kappa2(mp.pi/2), 25))
for key in ['pi_3','2pi_3','3pi_4','5pi_6']:
    a = ALPHA[key]
    print(key, 'kappa1=', mp.nstr(kappa1(a),20), ' kappa2=', mp.nstr(kappa2(a),16))
