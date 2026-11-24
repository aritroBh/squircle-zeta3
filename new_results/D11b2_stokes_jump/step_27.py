import json, time, math, sys
from mpmath import mp, mpf, mpc, pi, sqrt, nstr, exp, gamma, log
from mpmath.calculus.quadrature import GaussLegendre
from math import comb
tstart=time.time()
mp.dps = 120
sq2 = sqrt(2)
ipi = mpc(0,1)*pi
def pat(n): return [-1,1,1,-1][n%4]
GL = GaussLegendre(mp)
_nodes = {}
def glquad(f, a, b, n=96):
    if n not in _nodes:
        _nodes[n] = GL.calc_nodes(n, mp.dps+10)
    x, w = _nodes[n]
    s = mpc(0)
    h = (b-a)/2; m = (a+b)/2
    for xi, wi in zip(x, w):
        s += wi*f(m + h*xi)
    return h*s
def peval(c, t):
    v = c[0]
    for k in range(1,len(c)): v += c[k]*t**k
    return v
def t_of_w(w): return ipi*w/(1+w*w)
def tp_of_w(w): return ipi*(1-w*w)/(1+w*w)**2
pj = json.load(open('step_12_wpade.json'))
numR90 = [mpc(v) for v in pj['num']]; denR90 = [mpc(v) for v in pj['den']]
numR85 = [mpc(v) for v in pj['num85']]; denR85 = [mpc(v) for v in pj['den85']]
def PBWr90(w): return peval(numR90, w)/peval(denR90, w)
def PBWr85(w): return peval(numR85, w)/peval(denR85, w)
tp2 = 2/pi
NMAX=180
a2 = [mpf(0)] + [pat(n)*mpf(n)**mpf('-1.5')*tp2**n for n in range(1,NMAX+1)]
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
t0=time.time()
dm = d_series(a2)
numM90, denM90 = pade_coeffs(dm, 90, 90)
print("model pade built %.1fs"%(time.time()-t0)); sys.stdout.flush()
def PBWm90(w): return peval(numM90, w)/peval(denM90, w)
e34 = (-1+mpc(0,1))/sq2
def J_arc(PBW, eps, npts=96):
    f0 = math.sqrt(4*float(eps)/math.pi)
    phimax = min(1.30, 14*f0)
    def G(phi):
        w = mpc(math.cos(phi), math.sin(phi))
        return PBW(w)*tp_of_w(w)*mpc(0,1)*w
    def g(phi):
        return exp(-pi/(2*eps*math.cos(phi)))*(G(phi)-G(-phi))
    return glquad(g, mpf(0), mpf(repr(phimax)), npts)/(mpc(0,1)*eps)
def J_model_exact(eps):
    c0 = sq2/2*e34*2*pi*mpc(0,1)/gamma(mpf('1.5'))
    g = lambda y: exp(-y/eps)*mpc(0,1)*c0*sqrt(log(2*y/pi))
    v = glquad(g, pi/2, pi/2+30*eps, 96)
    return v/(mpc(0,1)*eps)
print("\neps   J_arc(real90)            J_arc(real85)            J_arc(model90)          J_model_exact         modelratio")
res = {}
for es in ['0.12','0.15','0.2','0.25','0.3']:
    eps = mpf(es)
    Jr90 = J_arc(PBWr90, eps)
    Jr85 = J_arc(PBWr85, eps)
    Jm = J_arc(PBWm90, eps)
    Je = J_model_exact(eps)
    res[es] = (Jr90, Jr85, Jm, Je)
    print("%s  %s  %s  %s  %s  %s" % (es, nstr(Jr90,11), nstr(Jr85,11), nstr(Jm,11), nstr(Je,11), nstr(Jm/Je,6)))
json.dump({es:[nstr(v,40) for v in res[es]] for es in res}, open('step_16_J.json','w'))
print("total %.1fs"%(time.time()-tstart))
