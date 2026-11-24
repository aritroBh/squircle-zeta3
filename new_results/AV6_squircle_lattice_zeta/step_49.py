code = r'''
from squircle import *
from methodC import *
from mpmath import mp, zeta, mpf, nstr, fabs, log, power, gamma, fac, pi
setdps(60)
s = 4
eps = mpf(1)/50
Z50 = ZA(50, 4)
diag = 8*zeta(s-1) - 4*zeta(s)*(1-power(2,-s*eps))
res = Z50 - diag
print("Z(50,4):", nstr(Z50, 40))
print("diag  :", nstr(diag, 40))
print("resid :", nstr(res, 12), "  resid/eps^4 =", nstr(res/eps**4, 15))
# predicted c_4 = (128/3) zeta(4) ln2^4 - 56 zeta(3) zeta(4)
u1 = gamma(s-1)*zeta(s-1)*zeta(s)*(1-power(2,1-s))
c4_pred = 4*zeta(s)*(s*log(2))**s/fac(s)*((-1)**s) - 8*s*u1
print("u1(4) =", nstr(u1, 20))
print("c4 pred =", nstr(c4_pred, 20))
# U_1 direct check of u1
def U1_direct(p, K=4000):
    tot = mpf(0)
    for k in range(2, K+1):
        inner = mpf(0)
        kf = mpf(k)
        for j in range(1, k):
            inner += log(1+power(j/kf, p))
        tot += power(kf, -s)*inner
    I0q = fsum([(-1)**(n+1)/(n*(n*p+1)) for n in range(1, 400)])
    tail = I0q*harmtail(K, s-1) - log(2)/2*harmtail(K, s)
    return tot + tail
for p_ in [40, 80]:
    u = U1_direct(p_)
    print(f"U1({p_})*p^3 = {nstr(u*p_**3, 15)}  (target {nstr(u1,12)})")
'''
exec(code)
with open("progress.log","a") as f: f.write("step_35: eps-structure: c_s formula, u_1(s) modular-Mellin exact form\n")
