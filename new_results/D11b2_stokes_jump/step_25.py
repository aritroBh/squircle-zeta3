import time, math, sys
from mpmath import mp, mpf, mpc, pi, sqrt, nstr, exp
from math import comb
mp.dps = 50
sq2 = sqrt(2)
def pat(n): return [-1,1,1,-1][n%4]
tp = 2/pi; ipi = mpc(0,1)*pi
NMAX=180
a2 = [mpf(0)] + [pat(n)*mpf(n)**mpf('-1.5')*tp**n for n in range(1,NMAX+1)]
def d_series(b):
    d = [mpc(0)]*(len(b)); d[0] = b[0]
    for k in range(1, len(b)):
        s = mpc(0)
        for m in range(0, k//2+1):
            n = k-2*m
            if n == 0: continue
            s += b[n] * ipi**n * (-1)**m * comb(k-m-1, m)
        d[k] = s
    return d
dm = d_series(a2)
def pade_coeffs(a, L, M):
    A = mp.matrix(M, M); rhs = mp.matrix(M, 1)
    for i in range(M):
        for j in range(M):
            k2 = L + i - j
            A[i, j] = a[k2] if k2 >= 0 else mpf(0)
        rhs[i] = -a[L + 1 + i]
    q = (A**-1)*rhs
    den = [mpf(1)] + [q[j] for j in range(M)]
    num = [sum(den[j]*a[n-j] for j in range(0, min(M, n)+1)) for n in range(L+1)]
    return num, den
numM, denM = pade_coeffs(dm, 90, 90)
def peval(c, t):
    v = c[0]
    for k in range(1,len(c)): v += c[k]*t**k
    return v
def PBWm(w): return peval(numM, w)/peval(denM, w)
def tp_of_w(w): return ipi*(1-w*w)/(1+w*w)**2
eps = mpf('0.2')
def G(phi):
    w = mpc(math.cos(phi), math.sin(phi))
    return PBWm(w)*tp_of_w(w)*mpc(0,1)*w
def g(phi):
    return exp(-pi/(2*eps*math.cos(phi)))*(G(phi)-G(-phi))
worst = []
for i in range(0, 200):
    phi = 1.5707963*i/199
    v = g(phi)
    worst.append((abs(v), phi))
worst.sort(reverse=True)
print("top |g| samples:", [(nstr(v,3), round(p,4)) for v,p in worst[:5]])
# scan near pi/2 closely
mx = 0; argmx=None
for i in range(1, 1000):
    phi = pi/2 - mpf(10)**(-i/50)
    v = abs(g(phi))
    if v>mx: mx=v; argmx=phi
print("max |g| near pi/2:", nstr(mx,3), "at pi/2-phi =", nstr(pi/2-argmx,3))
# also check |denM| minimum on the arc
mn = mpf(1e9); argmn=None
for i in range(0, 2000):
    phi = 1.5707963*i/1999
    w = mpc(math.cos(phi), math.sin(phi))
    dvv = abs(peval(denM, w))
    if dvv<mn: mn=dvv; argmn=phi
print("min |denM(e^{i phi})| =", nstr(mn,3), "at phi=", nstr(argmn,6))
