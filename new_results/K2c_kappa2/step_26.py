
import mpmath
from mpmath import mp, mpf, mpc, sqrt, log, exp, pi, gamma
import numpy as np
mp.dps = 140
u4 = exp(mpc(0,1)*mp.pi/4)

def gam(a): return mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4 - a/2))
def g1f(a): return -mp.mpf(3)/4 - 1j*mp.e**(1j*a)/(4*mp.sin(a))
def g2f(a):
    s = mp.sin(a); E = mp.e**(1j*a)
    return mp.mpf(29)/96 + E**2*(mp.mpf(1)/2 + 1/(32*s**2)) + (1j*E/s)*(E**2/4 - mp.mpf(3)/16)

def Jm(N, Z, m):
    a = abs(Z); s2 = mp.sqrt(2)
    f0 = lambda t: 1/s2 - (N/2)*(s2*a*t + 2*a*a)/(t**3 + s2*a*t*t + a*a*t)
    ts = mp.findroot(f0, mp.sqrt(N*a))
    xs = ts*u4
    Rp = ts/s2 + (N/2)*mp.log(1 + s2*a/ts + (a/ts)**2)   # Re phi at ts
    esc = exp(Rp)
    def ft(t):
        x = t*u4
        return exp(-x)*(1+Z/x)**(-N)*(x+Z)**(m+mpf(3)/2)/x * u4 * esc
    w = ts/6
    pts = [0, ts*mpf('0.3'), ts*mpf('0.6'), ts-w, ts-w/3, ts, ts+w/3, ts+w, ts*mpf('1.6'), ts*mpf('3'), ts*mpf('10'), mp.inf]
    return mp.quadgl(ft, pts)/esc

def Iplus(N, a):
    Z = 1j*a; tot = 0
    gms = [mp.mpf(1), g1f(a), g2f(a)]
    for m in range(3):
        bf = mp.gamma(m+mp.mpf('1.5'))*mp.gamma(N-m-mp.mpf('1.5'))/mp.gamma(N)
        tot += gms[m]*bf*Jm(N,Z,m)
    return (+1j/mp.pi)*gam(a)*tot

def beta1_formula(a):
    qq = mp.sqrt(1j*a)
    return 3*g1f(a)*qq/2 - qq**3/12 + 5*qq/4 + 1/(2*qq)
def beta2_formula(a):
    qq = mp.sqrt(1j*a); g1=g1f(a); g2=g2f(a)
    return (-g1*qq**4/8 + 21*g1*qq**2/8 + 9*g1/4 + 15*g2*qq**2/4
            + qq**6/288 - 5*qq**4/48 + 65*qq**2/96 + mp.mpf(15)/8)


a = mp.pi/2
A = mp.sqrt(2*a*mp.sin(a))
b1 = beta1_formula(a); b2f = beta2_formula(a)
print('N     sig      Re[(P-1-b1 s-b2 s^2)/s^3]   Im')
for N in [8000, 12000, 18000, 26000, 38000, 56000]:
    th = mp.sqrt(2*a*N)*(1+1j)
    P = -(2/A)*N*mp.e**th*Iplus(N, a)
    sig = 1/mp.sqrt(mpf(N))
    b3eff = (P-1-b1*sig-b2f*sig**2)/sig**3
    print('%6d %.5f  %s' % (N, float(sig), mp.nstr(b3eff, 18)))
