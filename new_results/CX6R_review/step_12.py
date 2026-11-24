import mpmath as mp, time
import cx6r_common as C
print("CX6R step_12: k=99 termwise |T_i| vs (12) summand; F(j,q;lam) vs (14)")
mp.mp.dps = 100
def betaL(L): return mp.log(1+mp.exp(-L))
def T_i(k, i, s, level=5):
    if i >= 1:
        phi = mp.arg(s+i+1); rot = mp.exp(-mp.j*phi); q = s+i+1
        def f(ell):
            L = ell*rot
            return mp.exp(-(s+1)*L)*(s*L)**(k-i)*betaL(L)**i/(1+mp.exp(-L))*rot
        cf = (-1)**i/(mp.factorial(i)*mp.factorial(k-i))
    else:
        phi = mp.arg(s+2); rot = mp.exp(-mp.j*phi); q = s+2
        def f(ell):
            L = ell*rot
            return mp.exp(-(s+2)*L)*(s*L)**k/(1+e**(-L) if False else mp.exp(-L))*rot
        cf = -1/mp.factorial(k)
    sc = 30/abs(q)
    segs = [[0, sc/3],[sc/3, sc],[sc, 3*sc],[3*sc, 8*sc],[8*sc, (k+40)/abs(q)+8*sc]]
    return cf*C.glint(f, segs, level)
def summand(k, i, s):
    rho = abs(s)
    if i == 0: return 2*rho**k/abs(s+2)**(k+1)
    return 2*(mp.mpf(2)**i/mp.factorial(i))*rho**(k-i)/abs(s+i+1)**(k-i+1)
k = 99
t = time.time()
for r in [mp.mpf('0'), mp.mpf('4'), mp.mpf('12')]:
    s = mp.mpf('0.5')+r*C.w4
    ratios = []
    for i in range(k+1):
        tv = abs(T_i(k, i, s))
        ratios.append((tv/summand(k, i, s), i))
    ratios.sort(reverse=True)
    top = ", ".join(f"i={i}: {mp.nstr(rt,4)}" for rt, i in ratios[:5])
    print(f"r={mp.nstr(r,4)}: worst |T_i|/summand ratios: {top}   ({time.time()-t:.0f}s)")
# F(j,q;lam) actual vs (14)
def Factual(j, q, lam):
    rs = mp.sqrt(j*q/(mp.sqrt(2)*lam)) if j > 0 else mp.mpf(1)
    pts = [0, rs/4, rs/2, rs, 2*rs, 4*rs, 10*rs, 10*rs+40/lam]
    segs = [[a,b] for a,b in zip(pts, pts[1:])]
    def f(r):
        s = mp.mpf('0.5')+r*C.w4
        return mp.exp(-lam*r)*abs(s)**j/abs(s+q)**(j+1)
    return C.glint(f, segs, 5)
def RHS14(j, q, lam):
    return 4/lam*mp.exp(lam/mp.mpf(2)+lam*q/mp.sqrt(2))*mp.exp(-mp.sqrt(3*lam*j*q/mp.sqrt(2)))
a = C.aa
worst = (mp.mpf(0), None)
t = time.time()
for i in range(0, 100):
    j = 99-i; q = i+1
    ratio = Factual(j, q, a)/RHS14(j, q, a)
    if ratio > worst[0]: worst = (ratio, (j, q))
print(f"F-bound k=99 lam=a: worst ratio {mp.nstr(worst[0],6)} at (j,q)={worst[1]}  ({time.time()-t:.0f}s)")
for kk in [99, 200, 400]:
    ratio = Factual(kk, 1, 3*a)/RHS14(kk, 1, 3*a)
    print(f"F-bound (j,q,lam)=({kk},1,3a): ratio {mp.nstr(ratio,6)}")
