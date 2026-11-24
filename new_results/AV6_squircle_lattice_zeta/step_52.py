code = r'''
from squircle import *
from methodC import *
from mpmath import mp, zeta, mpf, nstr, fabs, log, power, gamma, fac, matrix, lu_solve, pi
import json
setdps(65)
ps = [24, 32, 40, 48, 64, 96, 128]
epss = [mpf(1)/p_ for p_ in ps]
G = {}
G.update(json.load(open("epsgrid34.json")))
G.update(json.load(open("epsgrid56.json")))
def getZ(p_, s_): return mpf(G[f"{p_},{s_}"])

def lstsq(xs, ys, Q):
    # weighted LSQ for y = sum c_q x^q
    n = len(xs)
    A = matrix(n, Q+1)
    for i in range(n):
        xp = mpf(1)
        for q in range(Q+1):
            A[i,q] = xp; xp *= xs[i]
    # normal equations
    AtA = matrix(Q+1, Q+1); Atb = matrix(Q+1, 1)
    for a in range(Q+1):
        for b in range(Q+1):
            AtA[a,b] = sum(A[i,a]*A[i,b] for i in range(n))
        Atb[a,0] = sum(A[i,a]*ys[i] for i in range(n))
    c = lu_solve(AtA, Atb)
    resid = [ys[i] - sum(c[q,0]*xs[i]**q for q in range(Q+1)) for i in range(n)]
    return [c[q,0] for q in range(Q+1)], resid

for s_ in [3,4]:
    ys = [getZ(p_, s_) for p_ in ps]
    cs, resid = lstsq(epss, ys, 6)
    print(f"== s={s_}: max resid {nstr(max(abs(r) for r in resid),4)}")
    # predictions
    diag = lambda q: 4*zeta(s_)*((-1)**q)*(s_*log(2))**q/fac(q)
    print("c0:", nstr(cs[0], 30), " vs 8z(s-1):", nstr(8*zeta(s_-1), 30), " diff:", nstr(cs[0]-8*zeta(s_-1), 5))
    for q in range(1, s_):
        print(f"c{q}:", nstr(cs[q], 25), " vs diag:", nstr(diag(q), 25), " diff:", nstr(cs[q]-diag(q), 5))
    u1 = gamma(s_-1)*zeta(s_-1)*zeta(s_)*(1-power(2,1-s_))
    cs_pred = diag(s_) - 8*s_*u1
    print(f"c{s_}:", nstr(cs[s_], 25), " vs pred:", nstr(cs_pred, 25), " diff:", nstr(cs[s_]-cs_pred, 5))
    print(f"c{s_+1}:", nstr(cs[s_+1], 25))
    if s_+2 <= 6: print(f"c{s_+2}:", nstr(cs[s_+2], 25))
'''
exec(code)
