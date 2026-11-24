
import mpmath
from mpmath import mp, mpf, mpc, sqrt, log, exp, pi
import numpy as np
mp.dps = 70

def saddle_xs(N, Z):
    return (-Z + sqrt(Z**2 + 4*N*Z))/2

def Bbracket(m, N, Z):
    xs = saddle_xs(N, Z); w = xs + Z
    d2 = N*(1/xs**2 - 1/(xs+Z)**2)   # phi''(xs)
    ad2 = abs(d2)
    v = exp(-mpc(0,1)*mp.arg(d2)/2)   # steepest direction unit vector, v^2 d2 = |d2|
    phis = xs + N*log(1+Z/xs)
    def integrand(u):
        x = xs + u*v/sqrt(ad2)
        phi = x + N*log(1+Z/x)
        R = (xs/x) * ((x+Z)/w)**(m+mpf(3)/2)
        return R*exp(-(phi-phis))
    T = 16.0
    I = mp.quad(integrand, [-T, 0, T])
    return I/sqrt(2*pi)

def extract(m, Z, Ns):
    rows = []
    for N in Ns:
        sig = 1/sqrt(mpf(N))
        B = Bbracket(m, N, Z)
        rows.append((sig, (B-1)/sig, B))
    A = np.array([[sig**j for j in range(3)] for sig,_,_ in rows], dtype=complex)
    y = np.array([complex(v) for _,v,_ in rows], dtype=complex)
    sol, res, rk, sv = np.linalg.lstsq(A, y, rcond=None)
    return sol, rows

Z = mpc(0, mp.pi/2)
Ns = [1500, 2500, 4000, 6000, 9000, 14000]
for m in [0,1,2]:
    sol, rows = extract(m, Z, Ns)
    print('m=%d  b1=%.10f%+.10fi  b2=%.10f%+.10fi' % (m, sol[0].real, sol[0].imag, sol[1].real, sol[1].imag))
    for (sig, v, B) in rows:
        pred = sol[0]+sol[1]*sig+sol[2]*sig**2
        print('    sig=%.5f  B-1=%.8f%+.8fi  fitresid=%.2e' % (sig, float(mp.re(B-1)), float(mp.im(B-1)), abs(complex(v)-pred)))
print()
Zc = 1j*float(mp.pi)/2
print("K2b claimed numerics: b2: m=0: -0.625 ; m=1:", -1.75+3/(4*Zc), "; m=2:", -27/8+15/(4*Zc))
print("K2b claimed numerics: b1: m=0: 1/(2√Z)=%s  m=1: 3/(2√Z)=%s  m=2: 3/√Z=%s" % (1/(2*sqrt(Zc)), 3/(2*sqrt(Zc)), 3/sqrt(Zc)))
print("saddle5 symbolic  : b2: m=0:", -0.625+2835/(512*Zc), " m=1:", -1.75+699/(512*Zc), " m=2:", -27/8-285/(512*Zc))
