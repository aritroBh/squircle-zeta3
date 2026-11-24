
import mpmath
from mpmath import mp, mpf, mpc, sqrt, log, exp, pi
import numpy as np
mp.dps = 60

def saddle_xs(N, Z):
    return (-Z + sqrt(Z**2 + 4*N*Z))/2

def Bbracket(m, N, Z, Ucap=12.0):
    xs = saddle_xs(N, Z); w = xs + Z
    d2 = N*(1/xs**2 - 1/(xs+Z)**2)
    ad2 = abs(d2)
    v = exp(-mpc(0,1)*mp.arg(d2)/2)
    phis = xs + N*log(1+Z/xs)
    U = min(Ucap, 0.30*abs(xs)*sqrt(ad2))
    def integrand(u):
        x = xs + u*v/sqrt(ad2)
        phi = x + N*log(1+Z/x)
        R = (xs/x) * ((x+Z)/w)**(m+mpf(3)/2)
        return R*exp(-(phi-phis))
    I = mp.quadgl(integrand, [-U, 0, U])
    return I/sqrt(2*pi), U

def extract(m, Z, Ns):
    rows = []
    for N in Ns:
        sig = 1/sqrt(mpf(N))
        B, U = Bbracket(m, N, Z)
        rows.append((sig, (B-1)/sig, B, U))
    A = np.array([[sig**j for j in range(4)] for sig,_,_,_ in rows], dtype=complex)
    y = np.array([complex(v) for _,v,_,_ in rows], dtype=complex)
    sol, res, rk, sv = np.linalg.lstsq(A, y, rcond=None)
    return sol, rows

for Z, zname in [(mpc(0, mp.pi/2), 'i*pi/2'), (mpc(0,1), 'i')]:
    print('=== Z =', zname, '===')
    Ns = [160000, 260000, 420000, 680000, 1100000]
    for m in [0,1,2]:
        sol, rows = extract(m, Z, Ns)
        print('m=%d  b1=%.8f%+.8fi  b2=%.8f%+.8fi' % (m, sol[0].real, sol[0].imag, sol[1].real, sol[1].imag))
        for (sig, v, B, U) in rows:
            pred = sol[0]+sol[1]*sig+sol[2]*sig**2+sol[3]*sig**3
            print('    sig=%.6f U=%.1f  B-1=%.8f%+.8fi  fitresid=%.2e' % (sig, U, float(mp.re(B-1)), float(mp.im(B-1)), abs(complex(v)-pred)))
    Zc = complex(Z); sq = np.sqrt(Zc)
    print('  expect b1: m0 %s  m1 %s  m2 %s' % (1/(2*sq), 3/(2*sq), 3/sq))
    print('  K2b b2: m0 -0.625+0i; m1 %s; m2 %s' % (-1.75+3/(4*Zc), -27/8+15/(4*Zc)))
    print('  sad5 b2: m0 %s; m1 %s; m2 %s' % (-0.625+2835/(512*Zc), -1.75+699/(512*Zc), -27/8-285/(512*Zc)))
