
import mpmath
from mpmath import mp, mpf, mpc, sqrt, log, exp, pi
mp.dps = 60

def saddle_xs(N, Z):
    return (-Z + sqrt(Z**2 + 4*N*Z))/2

def Bbracket(m, N, Z):
    # exact loop integral for local moment m (Hankel power nu=m+1/2):
    # numerator J_m = int_0^inf e^{-phi(x)} x^{-1} (x+Z)^{m+3/2} dx, phi = x + N log(1+Z/x)
    xs = saddle_xs(N, Z); w = xs + Z
    d2 = N*(1/xs**2 - 1/(xs+Z)**2)   # phi''(xs)
    sig = 1/sqrt(mpf(N))
    width = 1/sqrt(abs(d2))
    a, b = xs - 14*width, xs + 14*width
    f = lambda x: exp(-(x + N*log(1+Z/x))) * (x+Z)**(m+mpf(3)/2) / x
    J = mp.quad(f, [a, xs, b])
    model = exp(-(xs + N*log(1+Z/xs))) * sqrt(2*pi/d2) * (xs+Z)**(m+mpf(3)/2)/xs
    return J/model

def extract(m, Z, Ns):
    rows = []
    for N in Ns:
        sig = 1/sqrt(mpf(N))
        B = Bbracket(m, N, Z)
        rows.append((sig, (B-1)/sig))
    # fit (B-1)/sig = b1 + b2*sig + b3*sig^2 by least squares
    A = [[sig**j for j in range(3)] for sig,_ in rows]
    y = [v for _,v in rows]
    # normal equations
    import itertools
    M = [[sum(A[k][i]*A[k][j] for k in range(len(rows))) for j in range(3)] for i in range(3)]
    rhs = [sum(A[k][i]*y[k] for k in range(len(rows))) for i in range(3)]
    # solve 3x3
    import numpy as np
    Mn = np.array([[float(mp.re(M[i][j])) for j in range(3)] for i in range(3)])
    rn = np.array([complex(rhs[i]) for i in range(3)], dtype=complex)
    # complex solve
    Mc = np.array([[complex(M[i][j]) for j in range(3)] for i in range(3)], dtype=complex)
    sol = np.linalg.solve(Mc, rn)
    return sol, rows

Z = mpc(0, mp.pi/2)
Ns = [1500, 2500, 4000, 6000, 9000, 14000]
for m in [0,1,2]:
    sol, rows = extract(m, Z, Ns)
    print('m=%d  b1=%.12f%+.12fi  b2=%.12f%+.12fi  b3=%.6f%+.6fi' % (
        m, sol[0].real, sol[0].imag, sol[1].real, sol[1].imag, sol[2].real, sol[2].imag))
    # residual check of fit
    for sig, v in rows:
        pred = sol[0]+sol[1]*sig+sol[2]*sig**2
        print('    sig=%.5f resid=%.2e' % (sig, abs(complex(v)-pred)))
print()
print("K2b claimed: b2^(0)=-5/8=-0.625; b2^(1)=-7/4+3/(4Z); b2^(2)=-27/8+15/(4Z)")
Zc = 1j*float(mp.pi)/2
print('  -7/4+3/(4Z) =', -1.75+3/(4*Zc))
print('  -27/8+15/(4Z) =', -27/8+15/(4*Zc))
print("saddle5 symbolic: m=0: -5/8+2835/(512Z) =", -0.625+2835/(512*Zc))
print("                  m=1: -7/4+699/(512Z) =", -1.75+699/(512*Zc))
print("                  m=2: -27/8-285/(512Z) =", -27/8-285/(512*Zc))
