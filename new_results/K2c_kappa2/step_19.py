
import mpmath
from mpmath import mp, mpf, mpc, sqrt, log, exp, pi, gamma
import numpy as np
mp.dps = 70

u4 = exp(mpc(0,1)*mp.pi/4)   # e^{i pi/4}

def gam(a): return mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4 - a/2))
def g1f(a): return -mp.mpf(3)/4 - 1j*mp.e**(1j*a)/(4*mp.sin(a))
def g2f(a):
    s = mp.sin(a); E = mp.e**(1j*a)
    return mp.mpf(29)/96 + E**2*(mp.mpf(1)/2 + 1/(32*s**2)) + (1j*E/s)*(E**2/4 - mp.mpf(3)/16)

def Jm(N, Z, m):
    a = abs(Z)
    # on diagonal x = t*u4: Re phi(t) = t/sqrt2 + (N/2) log(1 + sqrt2 a/t + a^2/t^2)
    s2 = mp.sqrt(2)
    def Rephi(t):
        return t/s2 + (N/2)*mp.log(1 + s2*a/t + a*a/t*t**0 + (a/t)**2)
    # find t* = argmin Rephi: solve d/dt = 0 : 1/s2 - (N/2)(s2 a t + 2 a^2)/(t^3 + s2 a t^2 + a^2 t) = 0
    f0 = lambda t: 1/s2 - (N/2)*(s2*a*t + 2*a*a)/(t**3 + s2*a*t*t + a*a*t)
    ts = mp.findroot(f0, mp.sqrt(N*a))
    def ft(t):
        x = t*u4
        return exp(-x)*(1+Z/x)**(-N)*(x+Z)**(m+mpf(3)/2)/x * u4
    pts = [0, ts*mpf('0.2'), ts*mpf('0.45'), ts*mpf('0.7'), ts*mpf('0.9'), ts,
           ts*mpf('1.12'), ts*mpf('1.3'), ts*mpf('1.7'), ts*mpf('3'), ts*mpf('8'), mp.inf]
    return mp.quadgl(ft, pts)

def Iplus(N, a):
    Z = 1j*a; tot = 0
    gms = [mp.mpf(1), g1f(a), g2f(a)]
    for m in range(3):
        bf = mp.gamma(m+mp.mpf('1.5'))*mp.gamma(N-m-mp.mpf('1.5'))/mp.gamma(N)
        tot += gms[m]*bf*Jm(N,Z,m)
    return (+1j/mp.pi)*gam(a)*tot

# sanity: J_0 vs saddle model
def saddle_xs(N, Z):
    return (-Z + sqrt(Z**2 + 4*N*Z))/2
for N in [1500, 10000]:
    Z = 1j*mp.pi/2; x1 = saddle_xs(N, Z)
    d2 = N*(1/x1**2 - 1/(x1+Z)**2)
    phis = x1 + N*log(1+Z/x1)
    model = exp(-phis)*sqrt(2*pi/d2)*(x1+Z)**(mpf(3)/2)/x1
    print('N=%d J0/model = %s' % (N, mp.nstr(Jm(N, Z, 0)/model, 12)))

def beta1_formula(a):
    qq = mp.sqrt(1j*a)
    return 3*g1f(a)*qq/2 - qq**3/12 + 5*qq/4 + 1/(2*qq)
def beta2_formula(a):
    qq = mp.sqrt(1j*a); g1=g1f(a); g2=g2f(a)
    return (-g1*qq**4/8 + 21*g1*qq**2/8 + 9*g1/4 + 15*g2*qq**2/4
            + qq**6/288 - 5*qq**4/48 + 65*qq**2/96 + mp.mpf(15)/8)

angles = [('pi/3', mp.pi/3), ('pi/2', mp.pi/2), ('2pi/3', 2*mp.pi/3), ('3pi/4', 3*mp.pi/4), ('5pi/6', 5*mp.pi/6)]
Ns = [1500, 2200, 3200, 4600, 6800, 10000]
print()
for name, a in angles:
    A = mp.sqrt(2*a*mp.sin(a))
    b1 = beta1_formula(a); b2f = beta2_formula(a)
    rows = []
    for N in Ns:
        th = mp.sqrt(2*a*N)*(1+1j)
        P = -(2/A)*N*mp.e**th*Iplus(N, a)
        sig = 1/mp.sqrt(mpf(N))
        rows.append((sig, P))
    Arows = [[sig**2, sig**3, sig**4] for sig, P in rows]
    yv = [complex(P - 1 - b1*sig) for sig, P in rows]
    sol, *_ = np.linalg.lstsq(np.array(Arows, dtype=complex), np.array(yv, dtype=complex), rcond=None)
    b2ext = mpc(sol[0].real, sol[0].imag)
    print('%s: |diff| = %s' % (name, mp.nstr(abs(b2ext-b2f), 3)))
    print('   extracted = %s' % mp.nstr(b2ext, 25))
    print('   formula   = %s' % mp.nstr(b2f, 25))
