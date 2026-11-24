
import mpmath
from mpmath import mp, mpf, mpc, sqrt, log, exp, pi, gamma, sin, cos
import numpy as np
mp.dps = 70

def gam(a): return mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4 - a/2))
def g1f(a): return -mp.mpf(3)/4 - 1j*mp.e**(1j*a)/(4*mp.sin(a))
def g2f(a):
    s = mp.sin(a); E = mp.e**(1j*a)
    return mp.mpf(29)/96 + E**2*(mp.mpf(1)/2 + 1/(32*s**2)) + (1j*E/s)*(E**2/4 - mp.mpf(3)/16)

def Jm(N, Z, m):
    xs = mp.sqrt(N*Z)
    f = lambda x: mp.e**(-x) * (1+Z/x)**(-N) * (x+Z)**(m+mp.mpf('1.5')) / x
    w = abs(xs)
    return mp.quad(f, [0, xs], maxdegree=10) + mp.quad(f, [xs, 2*xs, mp.inf], maxdegree=10)

def Iplus(N, a):
    Z = 1j*a; tot = 0
    gms = [mp.mpf(1), g1f(a), g2f(a)]
    for m in range(3):
        bf = mp.gamma(m+mp.mpf('1.5'))*mp.gamma(N-m-mp.mpf('1.5'))/mp.gamma(N)
        tot += gms[m]*bf*Jm(N,Z,m)
    return (+1j/mp.pi)*gam(a)*tot

def beta1_formula(a):
    Z = 1j*a; qq = mp.sqrt(Z)
    return 3*g1f(a)*qq/2 - qq**3/12 + 5*qq/4 + 1/(2*qq)

def beta2_formula(a):
    Z = 1j*a; qq = mp.sqrt(Z); g1=g1f(a); g2=g2f(a)
    return (-g1*qq**4/8 + 21*g1*qq**2/8 + 9*g1/4 + 15*g2*qq**2/4
            + qq**6/288 - 5*qq**4/48 + 65*qq**2/96 + mp.mpf(15)/8)

angles = [('pi/3', mp.pi/3), ('pi/2', mp.pi/2), ('2pi/3', 2*mp.pi/3), ('3pi/4', 3*mp.pi/4), ('5pi/6', 5*mp.pi/6)]
Ns = [2500, 3600, 5000, 7000, 10000]
for name, a in angles:
    A = mp.sqrt(2*a*mp.sin(a))
    b1 = beta1_formula(a)
    rows = []
    for N in Ns:
        th = mp.sqrt(2*a*N)*(1+1j)
        P = -(2/A)*N*mp.e**th*Iplus(N, a)
        sig = 1/mp.sqrt(mpf(N))
        rows.append((sig, P))
    # fit P = 1 + b1 sig + b2 sig^2 + b3 sig^3 + b4 sig^4 (b1 fixed to formula)
    Arows, yv = [], []
    for sig, P in rows:
        Arows.append([sig**2, sig**3, sig**4])
        yv.append(complex(P - 1 - b1*sig))
    Am = np.array(Arows, dtype=complex); ym = np.array(yv, dtype=complex)
    sol, *_ = np.linalg.lstsq(Am, ym, rcond=None)
    b2ext = mpc(sol[0].real, sol[0].imag)
    b2f = beta2_formula(a)
    print('%s: beta2 extracted = %s' % (name, mp.nstr(b2ext, 20)))
    print('%s: beta2 formula   = %s   |diff| = %s' % (name, mp.nstr(b2f, 20), mp.nstr(abs(b2ext-b2f), 3)))
    # fit quality
    for (sig, P), ar, yy in zip(rows, Arows, yv):
        pred = sol[0]*ar[0]+sol[1]*ar[1]+sol[2]*ar[2]
    print()
