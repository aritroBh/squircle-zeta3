
import mpmath
from mpmath import mp, mpf, mpc, sqrt, log, exp, pi, gamma
mp.dps = 120
u4 = exp(mpc(0,1)*mp.pi/4)
s2 = mp.sqrt(2)

def gam(a): return mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4 - a/2))
def g1f(a): return -mp.mpf(3)/4 - 1j*mp.e**(1j*a)/(4*mp.sin(a))
def g2f(a):
    s = mp.sin(a); E = mp.e**(1j*a)
    return mp.mpf(29)/96 + E**2*(mp.mpf(1)/2 + 1/(32*s**2)) + (1j*E/s)*(E**2/4 - mp.mpf(3)/16)

def Jm(N, Z, m, p):
    a = abs(Z); sq2 = mp.sqrt(2)
    f0 = lambda t: p/sq2 - (N/2)*(sq2*a*t + 2*a*a)/(t**3 + sq2*a*t*t + a*a*t)
    ts = mp.findroot(f0, mp.sqrt(N*a/p))
    Rp = ts*p/sq2 + (N/2)*mp.log(1 + sq2*a/ts + (a/ts)**2)
    esc = exp(Rp)
    def ft(t):
        x = t*u4
        return exp(-p*x)*(1+Z/x)**(-N)*(x+Z)**(m+mpf(3)/2)/x * u4 * esc
    w = ts/6
    pts = [0, ts*mpf('0.3'), ts*mpf('0.6'), ts-w, ts-w/3, ts, ts+w/3, ts+w, ts*mpf('1.6'), ts*mpf('3'), ts*mpf('10'), mp.inf]
    return mp.quadgl(ft, pts)/esc

def Iplus_m(N, a, p, m):
    Z = 1j*a
    bf = mp.gamma(m+mp.mpf('1.5'))*mp.gamma(N-m-mp.mpf('1.5'))/mp.gamma(N)
    gms = [mp.mpf(1), g1f(a), g2f(a)]
    return (+1j/mp.pi)*gam(a)*gms[m]*bf*Jm(N,Z,m,p)

def Iplus(N, a, p):
    return sum(Iplus_m(N,a,p,m) for m in range(3))


import numpy as np
# d1 from symbolic formula (step_32): d1 = (q^2*(-9 sqrt2 g1 + 18 g1 - q^2 + 2 sqrt2 q^2 - 15 sqrt2 + 15) - 3 sqrt2 + 6)/(12 q), tot = -(1+d1 sig)
def d1_formula(Z, g1):
    q = mp.sqrt(Z); sq2 = mp.sqrt(2)
    return -( (q**2*(-9*sq2*g1 + 18*g1 - q**2 + 2*sq2*q**2 - 15*sq2 + 15) - 3*sq2 + 6)/(12*q) )
    # minus because totn = -(1+d1 sig): C_img = -e^{Z/2}(1 + d1 sig) with d1 = -d1_printed

for a in [mp.pi/3, mp.pi/2, 2*mp.pi/3, 3*mp.pi/4, 5*mp.pi/6]:
    Z = 1j*a; g1 = g1f(a)
    rows = []
    for N in [1500, 3000, 6000, 12000]:
        sig = 1/mp.sqrt(mpf(N))
        th = mp.sqrt(2*a*N)
        I1 = Iplus(N, a, 1); I2 = Iplus(N, a, 2)
        C = -2*I2/I1*exp((s2-1)*th*(1+1j))
        rows.append((sig, C/(-exp(Z/2))))
    Arows = [[s, s*s] for s,_ in rows]
    yv = [complex(v-1) for _,v in rows]
    sol, *_ = np.linalg.lstsq(np.array(Arows,dtype=complex), np.array(yv,dtype=complex), rcond=None)
    d1ext = mpc(sol[0].real, sol[0].imag)
    d1f = d1_formula(Z, g1)
    print('alpha=%s' % mp.nstr(a,5))
    print('  d1 extracted = %s   formula = %s   |diff|=%s' % (mp.nstr(d1ext,16), mp.nstr(d1f,16), mp.nstr(abs(d1ext-d1f),3)))
