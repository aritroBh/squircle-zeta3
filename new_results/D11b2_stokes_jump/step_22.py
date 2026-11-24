import json, time, math
from mpmath import mp, mpf, mpc, pi, factorial, sqrt, nstr, quad, exp, polylog
from math import comb
t0=time.time()
mp.dps = 100
sq2 = sqrt(2)
def pat(n): return [-1,1,1,-1][n%4]
tp = 2/pi
ipi = mpc(0,1)*pi
NMAX=180
a2 = [mpf(0)] + [pat(n)*mpf(n)**mpf('-1.5')*tp**n for n in range(1,NMAX+1)]
e34 = (-1+mpc(0,1))/sq2   # e^{3i pi/4}
e34m = (-1-mpc(0,1))/sq2  # e^{-3i pi/4}
s32 = mpf('1.5')
def Bm(t):
    return sq2/2*( e34*polylog(s32, -2*mpc(0,1)*t/pi) + e34m*polylog(s32, 2*mpc(0,1)*t/pi) )
print("Bm(0.3) vs series:", nstr(abs(Bm(mpf('0.3')) - sum(a2[n]*mpf('0.3')**n for n in range(NMAX+1))),3))
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
def BmW(w): return Bm(t_of_w(w))
print("\nw-Pade[90/90] model vs exact on arcs:")
for phi in [0.0, 0.3, 0.6, 0.9, 1.2]:
    w = mpc(math.cos(phi), math.sin(phi))
    pv = PBWm(w); ev = BmW(w)
    print("  phi=%.1f  |pade-exact| = %s" % (phi, nstr(abs(pv-ev),3)))

# ---- exact model jump: J = (1/(i eps)) int_{pi/2}^inf e^{-y/eps} i [Bm(iy-0) - Bm(iy+0)] dy
# left side Re t<0: t = iy - sig, sig>0 ; z = -2it/pi = (2y+2i sig)/pi = X + i0
def J_model_exact(eps):
    def g(y):
        X = 2*y/pi
        zl = mpc(X, 1e-8)   # left (Re t<0) -> z = X + i0
        zr = mpc(X, -1e-8)
        disc = sq2/2*e34*(polylog(s32, zl) - polylog(s32, zr))
        return exp(-y/eps)*mpc(0,1)*disc
    v = quad(g, [pi/2, pi/2+2*eps, pi/2+10*eps, pi/2+40*eps])
    return v/(mpc(0,1)*eps)
def J_arc(PBW, eps, phimax=pi/2):
    def G(phi):
        w = mpc(math.cos(phi), math.sin(phi))
        return PBW(w)*tp_of_w(w)*mpc(0,1)*w
    def g(phi):
        return exp(-pi/(2*eps*math.cos(phi)))*(G(phi)-G(-phi))
    f0 = math.sqrt(4*eps/float(pi))
    v = quad(g, [0, f0, 3*f0, 8*f0, phimax])
    return v/(mpc(0,1)*eps)
print("\nmodel J: exact Hankel vs arc-exact vs arc-Pade:")
for es in ['0.15','0.2','0.25','0.3']:
    eps = mpf(es)
    Je = J_model_exact(eps)
    Ja = J_arc(BmW, eps)
    Jp = J_arc(PBWm, eps)
    print(f"  eps={es}: exact={nstr(Je,10)}")
    print(f"           arcex={nstr(Ja,10)}  rel diff={nstr(abs(Ja-Je)/abs(Je),2)}")
    print(f"           pade ={nstr(Jp,10)}  rel diff={nstr(abs(Jp-Je)/abs(Je),2)}")
print("(%.1fs)"%(time.time()-t0))
