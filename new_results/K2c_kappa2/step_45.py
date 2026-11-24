work='<project>/research/08_theory/K2c_kappa2/'
code = r'''
import sys; sys.path.insert(0,'<project>/research/08_theory/K2b_corner_phase/')
from k2utils import load
import mpmath
from mpmath import mp, mpf, mpc, sqrt, log, exp, pi
mp.dps = 80
s2 = mp.sqrt(2)

def g1f(a): return -mp.mpf(3)/4 - 1j*mp.e**(1j*a)/(4*mp.sin(a))
def kap1(a):
    C = 1/(2*mp.sqrt(a)) + a**mpf('1.5')/12 + 3*mp.sqrt(a)/8/mp.tan(a)
    return (C+mp.sqrt(a)/2)/s2 + 1j*(C-mp.sqrt(a)/2)/s2
def kap2(a):
    return (mpf(3)/4 + a*a/24 + 3*a/(16*mp.tan(a))
            + 1j*(a**3/288 - 11*a/384 + (a*a/32 + mpf(9)/16)/mp.tan(a) - 15*a/(128*mp.sin(a)**2)))
def beta1(a):
    qq = mp.sqrt(1j*a); g1 = g1f(a)
    return 3*g1*qq/2 - qq**3/12 + 5*qq/4 + 1/(2*qq)
def beta2(a): return mp.conj(kap2(a))
def d1f(a):
    Z = 1j*a; q = mp.sqrt(Z); g1 = g1f(a); sq2 = s2
    return -( (q**2*(-9*sq2*g1 + 18*g1 - q**2 + 2*sq2*q**2 - 15*sq2 + 15) - 3*sq2 + 6)/(12*q) )

angles = [('pi_3', mp.pi/3), ('pi_2', mp.pi/2), ('2pi_3', 2*mp.pi/3), ('3pi_4', 3*mp.pi/4), ('5pi_6', 5*mp.pi/6)]
print('image amplitude check: -e^{i a/2} at pi/2 =', mp.nstr(-exp(1j*(mp.pi/2)/2),12))
print()
for key, a in angles:
    d = load(key)
    A = mp.sqrt(2*a*mp.sin(a))
    k1, k2, b1, b2, d1 = kap1(a), kap2(a), beta1(a), beta2(a), d1f(a)
    print('=== alpha = %s ===  kappa2 = %s' % (key, mp.nstr(k2, 16)))
    print('  k      rho2*N^1.5        rho2-img*N^1.5   (img model: -e^{Za/2} e^{-(s2-1)th(1+i)}(1+d1 sig), d1=%s)' % mp.nstr(d1,8))
    for k in [100, 125, 150, 175, 200, 225, 250, 275, 300]:
        N = mpf(k+1); th = mp.sqrt(2*a*N); sig = 1/mp.sqrt(N)
        W = -mp.re(d[k])*N*exp(th)/A
        c, s = mp.cos(th), mp.sin(th)
        def ReE(v): return mp.re(v)*c - mp.im(v)*s   # Re[e^{i th} v]
        pred = ReE(1) + ReE(k1)*sig + ReE(k2)*sig**2
        rho2 = W - pred
        # image in W: Re[e^{-i th} P R], P = 1+b1 sig+b2 sig^2, R = -e^{Z/2} e^{-(s2-1)th(1+i)}(1+d1 sig)
        P = 1 + b1*sig + b2*sig**2
        R = -exp(1j*a/2)*exp(-(s2-1)*th*(1+1j))*(1+d1*sig)
        Wimg = mp.re(exp(-1j*th)*P*R)
        print('  %3d  %+14s  %+14s' % (k, mp.nstr(rho2*N**mpf('1.5'), 12), mp.nstr((rho2-Wimg)*N**mpf('1.5'), 12)))
    print()
'''
open(work+'step_34.py','w').write(code)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python', work+'step_34.py'],capture_output=True,text=True,timeout=1200)
print(r.stdout); print(r.stderr[-500:])
