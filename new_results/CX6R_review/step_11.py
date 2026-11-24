import mpmath as mp, time
import cx6r_common as C
print("CX6R step_11: rotated finite-sum form of A_k (sec.5), validity + termwise bounds")
mp.mp.dps = 100
def betaL(L): return mp.log(1+mp.exp(-L))
def T_i(k, i, s, level=5):
    # i>=1:  (-1)^i/(i!(k-i)!) int_0^inf e^{-(s+1)L} (sL)^{k-i} beta(L)^i/(1+e^{-L}) dL
    #        rotated: L = ell * exp(-i phi_i), phi_i = arg(s+i+1)
    # i=0:  -(1/k!) int e^{-(s+2)L} (sL)^k/(1+e^{-L}) dL, rotated by arg(s+2)
    if i >= 1:
        phi = mp.arg(s+i+1); rot = mp.exp(-mp.j*phi)
        q = s+i+1
        def f(ell):
            L = ell*rot
            return mp.exp(-(s+1)*L)*(s*L)**(k-i)*betaL(L)**i/(1+mp.exp(-L))*rot
        cf = (-1)**i/(mp.factorial(i)*mp.factorial(k-i))
    else:
        phi = mp.arg(s+2); rot = mp.exp(-mp.j*phi)
        q = s+2
        def f(ell):
            L = ell*rot
            return mp.exp(-(s+2)*L)*(s*L)**k/(1+mp.exp(-L))*rot
        cf = -1/mp.factorial(k)
    sc = 30/abs(q)
    segs = [[0, sc/3],[sc/3, sc],[sc, 3*sc],[3*sc, 8*sc],[8*sc, (k+40)/abs(q)+8*sc]]
    return cf*C.glint(f, segs, level)
def summand(k, i, s):
    rho = abs(s)
    if i == 0: return 2*rho**k/abs(s+2)**(k+1)
    return 2*(mp.mpf(2)**i/mp.factorial(i))*rho**(k-i)/abs(s+i+1)**(k-i+1)
def Ak_sum(k, s, level=5):
    return sum(T_i(k, i, s, level) for i in range(k+1))
# 1) validity of expansion+rotation at moderate k (direct A_k available)
for k in [5, 10, 20, 30]:
    for r in [mp.mpf('0'), mp.mpf('2'), mp.mpf('8')]:
        s = mp.mpf('0.5')+r*C.w4
        direct = C.Ak(k, s, level=6)
        viasum = Ak_sum(k, s)
        print(f"k={k:2d} r={mp.nstr(r,3)}: |direct - Sum T_i| = {mp.nstr(abs(direct-viasum),5)}")
