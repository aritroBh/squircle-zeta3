import json, time, math
from mpmath import mp, mpf, mpc, pi, factorial, sqrt, nstr, quad, exp, gamma, log
from math import comb
t0=time.time()
mp.dps = 70
sq2 = sqrt(2)
def pat(n): return [-1,1,1,-1][n%4]
tp = 2/pi
ipi = mpc(0,1)*pi
NMAX=180
a2 = [mpf(0)] + [pat(n)*mpf(n)**mpf('-1.5')*tp**n for n in range(1,NMAX+1)]
e34 = (-1+mpc(0,1))/sq2
def t_of_w(w): return ipi*w/(1+w*w)
def tp_of_w(w): return ipi*(1-w*w)/(1+w*w)**2
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
# exact model jump via analytic discontinuity: disc Li_{3/2}(X) = 2 pi i (ln X)^{1/2}/Gamma(3/2)
def J_model_exact(eps):
    c0 = sq2/2*e34*2*pi*mpc(0,1)/gamma(mpf('1.5'))
    def g(y):
        X = 2*y/pi
        return exp(-y/eps)*mpc(0,1)*c0*sqrt(log(X))
    v = quad(g, [pi/2, pi/2+3*eps, pi/2+20*eps, mp.inf])
    return v/(mpc(0,1)*eps)
def J_arc(PBW, eps):
    def G(phi):
        w = mpc(math.cos(phi), math.sin(phi))
        return PBW(w)*tp_of_w(w)*mpc(0,1)*w
    def g(phi):
        return exp(-pi/(2*eps*math.cos(phi)))*(G(phi)-G(-phi))
    f0 = math.sqrt(4*float(eps)/math.pi)
    v = quad(g, [0, f0, 3*f0, 9*f0, pi/2])
    return v/(mpc(0,1)*eps)
print("model J: exact-Hankel vs arc-Pade   (prediction: 2 sqrt(pi) e^{-3i pi/4} eps^{1/2} e^{-pi/(2 eps)})")
for es in ['0.12','0.15','0.2','0.25','0.3']:
    eps = mpf(es)
    Je = J_model_exact(eps)
    Jp = J_arc(PBWm, eps)
    Jpred = 2*sqrt(pi)*(-1-mpc(0,1))/sq2 * sqrt(eps)*exp(-pi/(2*eps))
    print(f"  eps={es}: exact={nstr(Je,12)}")
    print(f"       arcPade={nstr(Jp,12)}  rel diff={nstr(abs(Jp-Je)/abs(Je),2)}   exact/pred={nstr(Je/Jpred,8)}")
print("(%.1fs)"%(time.time()-t0))
