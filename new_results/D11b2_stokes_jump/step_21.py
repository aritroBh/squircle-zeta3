import json, time
from mpmath import mp, mpf, mpc, pi, factorial, sqrt, nstr, quad, exp, polylog
from math import comb
t0=time.time()
mp.dps = 100
sq2 = sqrt(2)
def pat(n): return [-1,1,1,-1][n%4]
tp = 2/pi
ipi = mpc(0,1)*pi
# model Borel coeffs (leading D11 asymptotic, exact): a_n = pat(n) n^{-3/2} (2/pi)^n
NMAX=180
a2 = [mpf(0)] + [pat(n)*mpf(n)**mpf('-1.5')*tp**n for n in range(1,NMAX+1)]
# exact model function: Bm(t) = (sq2/2)[e^{3i pi/4} Li_{3/2}(-2it/pi) + e^{-3i pi/4} Li_{3/2}(2it/pi)]
cp = sq2/2*mpc(1,1)/sq2  # e^{3i pi/4} = (1+i)/sqrt2
e34 = (1+mpc(0,1))/sq2   # e^{3i pi/4}
e34m = (1-mpc(0,1))/sq2  # e^{-3i pi/4}
def Bm(t):
    return sq2/2*( e34*polylog(mpf('1.5'), -2*mpc(0,1)*t/pi) + e34m*polylog(mpf('1.5'), 2*mpc(0,1)*t/pi) )
# series check
print("Bm series vs a2: ", nstr(abs(Bm(mpf('0.3')) - sum(a2[n]*mpf('0.3')**n for n in range(NMAX+1))),3))

def t_of_w(w): return ipi*w/(1+w*w)
def tp_of_w(w): return ipi*(1-w*w)/(1+w*w)**2
def d_series(b):
    d = [mpc(0)]*(len(b))
    d[0] = b[0]
    for k in range(1, len(b)):
        s = mpc(0)
        for m in range(0, k//2+1):
            n = k-2*m
            if n == 0: continue
            s += b[n] * ipi**n * (-1)**m * comb(k-m-1, m)
        d[k] = s
    return d
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
def peval(c, t):
    v = c[0]
    for k in range(1,len(c)): v += c[k]*t**k
    return v

dm = d_series(a2)
numM, denM = pade_coeffs(dm, 90, 90)
def PBWm(w): return peval(numM, w)/peval(denM, w)
def BmW(w): return Bm(t_of_w(w))   # exact mapped model (one branch = inside disc)

# check w-Pade of model against exact mapped model inside disc and on arcs
print("\nw-Pade[90/90] of model vs exact Bm(t(w)):")
for phi in [0, 0.3, 0.6, 0.9, 1.2]:
    w = mpc(*__import__('math').polar(1, phi)) if False else mpc(__import__('math').cos(phi), __import__('math').sin(phi))
    pv = PBWm(w); ev = BmW(w)
    print("  phi=%.1f  pade=%s  exact=%s  diff=%s" % (phi, nstr(pv,8), nstr(ev,8), nstr(abs(pv-ev),2)))
print("(%.1fs)"%(time.time()-t0))
