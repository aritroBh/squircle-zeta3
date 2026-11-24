code = r'''import mpmath as mp, time
import cx6r_common as C
print("CX6R step_12b: k=99 termwise + F-bound, HIGH GL levels with convergence check")
mp.mp.dps = 100
def betaL(L): return mp.log(1+mp.exp(-L))
def T_i(k, i, s, level):
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
            return mp.exp(-(s+2)*L)*(s*L)**k/(1+mp.exp(-L))*rot
        cf = -1/mp.factorial(k)
    sc = 30/abs(q)
    segs = [[0, sc/3],[sc/3, sc],[sc, 3*sc],[3*sc, 8*sc],[8*sc, (k+40)/abs(q)+8*sc]]
    return cf*C.glint(f, segs, level)
def summand(k, i, s):
    rho = abs(s)
    if i == 0: return 2*rho**k/abs(s+2)**(k+1)
    return 2*(mp.mpf(2)**i/mp.factorial(i))*rho**(k-i)/abs(s+i+1)**(k-i+1)
k = 99
# convergence check of T_0 at s=1/2 (exact value available)
s = mp.mpf('0.5')
exact_T0 = -s**k*mp.nsum(lambda n: (-1)**n/(s+2+n)**(k+1), [0, mp.inf])
for lev in (5,6,7,8):
    v = T_i(k, 0, s, lev)
    print(f"T_0 level {lev}: |err vs exact| = {mp.nstr(abs(v-exact_T0),5)}")
print("exact |T_0|/summand =", mp.nstr(abs(exact_T0)/summand(k,0,s), 10))
t=time.time()
for r in [mp.mpf('0'), mp.mpf('4'), mp.mpf('12')]:
    s = mp.mpf('0.5')+r*C.w4
    ratios = []
    for i in range(k+1):
        tv = abs(T_i(k, i, s, 7))
        ratios.append((tv/summand(k, i, s), i))
    ratios.sort(reverse=True)
    top = ", ".join(f"i={i}: {mp.nstr(rt,4)}" for rt, i in ratios[:5])
    print(f"r={mp.nstr(r,4)}: worst |T_i|/summand: {top}  ({time.time()-t:.0f}s)")
# F bound with convergence check
def Factual(j, q, lam, level):
    rs = mp.sqrt(j*q/(mp.sqrt(2)*lam)) if j > 0 else mp.mpf(1)
    pts = [0, rs/4, rs/2, rs, 2*rs, 4*rs, 10*rs, 10*rs+40/lam]
    segs = [[a,b] for a,b in zip(pts, pts[1:])]
    def f(r):
        s = mp.mpf('0.5')+r*C.w4
        return mp.exp(-lam*r)*abs(s)**j/abs(s+q)**(j+1)
    return C.glint(f, segs, level)
def RHS14(j, q, lam):
    return 4/lam*mp.exp(lam/mp.mpf(2)+lam*q/mp.sqrt(2))*mp.exp(-mp.sqrt(3*lam*j*q/mp.sqrt(2)))
a = C.aa
v6 = Factual(99, 1, a, 6); v7 = Factual(99, 1, a, 7); v8 = Factual(99, 1, a, 8)
print("F(99,1,a) levels 6,7,8:", mp.nstr(v6,12), mp.nstr(v7,12), mp.nstr(v8,12))
t=time.time()
worst = (mp.mpf(0), None)
for i in range(0, 100):
    j = 99-i; q = i+1
    ratio = Factual(j, q, a, 7)/RHS14(j, q, a)
    if ratio > worst[0]: worst = (ratio, (j, q))
print(f"F-bound k=99 lam=a: worst ratio {mp.nstr(worst[0],6)} at (j,q)={worst[1]}  ({time.time()-t:.0f}s)")
for kk in [99, 200, 400]:
    r7 = Factual(kk, 1, 3*a, 7); r8 = Factual(kk, 1, 3*a, 8)
    print(f"F({kk},1,3a): ratio(lev7) {mp.nstr(r7/RHS14(kk,1,3*a),6)}  ratio(lev8) {mp.nstr(r8/RHS14(kk,1,3*a),6)}")
'''
open('step_12b.py','w').write(code)
exec(compile(code,'step_12b.py','exec'))
