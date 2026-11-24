
import mpmath
from mpmath import mp, mpf, mpc, sqrt, log, exp, pi, gamma
import numpy as np
mp.dps = 60

def gam(a): return mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4 - a/2))
def g1f(a): return -mp.mpf(3)/4 - 1j*mp.e**(1j*a)/(4*mp.sin(a))
def g2f(a):
    s = mp.sin(a); E = mp.e**(1j*a)
    return mp.mpf(29)/96 + E**2*(mp.mpf(1)/2 + 1/(32*s**2)) + (1j*E/s)*(E**2/4 - mp.mpf(3)/16)

def saddle_xs(N, Z):
    return (-Z + sqrt(Z**2 + 4*N*Z))/2

def Jm(N, Z, m):
    x1 = saddle_xs(N, Z)
    def f(x):
        return exp(-x)*(1+Z/x)**(-N)*(x+Z)**(m+mpf(3)/2)/x
    # piece 1: 0 -> x1  (x = t*x1, t in [0,1])
    I1 = x1*mp.quadgl(lambda t: f(t*x1), [0, mpf('0.2'), mpf('0.5'), mpf('0.8'), 1])
    # piece 2: x1 -> x1+T*v (steepest dir), T so that Re(Dphi) ~ 60
    d2 = N*(1/x1**2 - 1/(x1+Z)**2); ad2 = abs(d2)
    v = exp(-mpc(0,1)*mp.arg(d2)/2)
    T = 11.0/sqrt(ad2)
    x2 = x1 + T*v
    I2 = (x2-x1)*mp.quadgl(lambda t: f(x1+t*(x2-x1)), [0, mpf('0.3'), mpf('0.7'), 1])
    # piece 3: x2 -> x2 + inf along real direction
    I3 = mp.quadgl(lambda s: f(x2+s), [0, mpf('5'), mpf('20'), mp.inf])
    return I1+I2+I3

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
Ns = [1500, 2200, 3200, 4600, 6800, 10000]
out = {}
for name, a in angles:
    A = mp.sqrt(2*a*mp.sin(a))
    b1 = beta1_formula(a)
    rows = []
    for N in Ns:
        th = mp.sqrt(2*a*N)*(1+1j)
        P = -(2/A)*N*mp.e**th*Iplus(N, a)
        sig = 1/mp.sqrt(mpf(N))
        rows.append((sig, P))
    Arows = [[sig**2, sig**3, sig**4, sig**5] for sig, P in rows]
    yv = [complex(P - 1 - b1*sig) for sig, P in rows]
    sol, *_ = np.linalg.lstsq(np.array(Arows, dtype=complex), np.array(yv, dtype=complex), rcond=None)
    b2ext = mpc(sol[0].real, sol[0].imag)
    b2f = beta2_formula(a)
    d = abs(b2ext-b2f)
    out[name] = (b2ext, b2f, d)
    print('%s: beta2 extracted = %s' % (name, mp.nstr(b2ext, 22)))
    print('%s: beta2 formula   = %s   |diff| = %s' % (name, mp.nstr(b2f, 22), mp.nstr(d, 3)))
    # residual of fit at each N
    for (sig, P), ar, yy in zip(rows, Arows, yv):
        pred = sum(sol[j]*ar[j] for j in range(4))
    print()
