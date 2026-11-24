import time, math, sys
from mpmath import mp, mpf, mpc, pi, sqrt, nstr, quad, exp, gamma, log
from math import comb
mp.dps = 50
sq2 = sqrt(2)
def pat(n): return [-1,1,1,-1][n%4]
tp = 2/pi; ipi = mpc(0,1)*pi
NMAX=180
t0=time.time()
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
print("d_series %.2fs"%(time.time()-t0)); sys.stdout.flush()
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
t0=time.time()
numM, denM = pade_coeffs(dm, 90, 90)
print("pade %.1fs"%(time.time()-t0)); sys.stdout.flush()
def peval(c, t):
    v = c[0]
    for k in range(1,len(c)): v += c[k]*t**k
    return v
def PBWm(w): return peval(numM, w)/peval(denM, w)
e34 = (-1+mpc(0,1))/sq2
def J_model_exact(eps):
    c0 = sq2/2*e34*2*pi*mpc(0,1)/gamma(mpf('1.5'))
    def g(y):
        X = 2*y/pi
        return exp(-y/eps)*mpc(0,1)*c0*sqrt(log(X))
    v = quad(g, [pi/2, pi/2+3*eps, pi/2+20*eps])
    return v/(mpc(0,1)*eps)
t0=time.time()
Je = J_model_exact(mpf('0.2'))
print("Jexact(0.2) =", nstr(Je,12), "  %.1fs"%(time.time()-t0)); sys.stdout.flush()
def t_of_w(w): return ipi*w/(1+w*w)
def tp_of_w(w): return ipi*(1-w*w)/(1+w*w)**2
def J_arc(PBW, eps):
    def G(phi):
        w = mpc(math.cos(phi), math.sin(phi))
        return PBW(w)*tp_of_w(w)*mpc(0,1)*w
    def g(phi):
        return exp(-pi/(2*eps*math.cos(phi)))*(G(phi)-G(-phi))
    f0 = math.sqrt(4*float(eps)/math.pi)
    v = quad(g, [0, f0, 3*f0, 9*f0, pi/2])
    return v/(mpc(0,1)*eps)
t0=time.time()
Jp = J_arc(PBWm, mpf('0.2'))
print("JarcPade(0.2) =", nstr(Jp,12), "  %.1fs"%(time.time()-t0)); sys.stdout.flush()
print("ratio:", nstr(Jp/Je, 10))
