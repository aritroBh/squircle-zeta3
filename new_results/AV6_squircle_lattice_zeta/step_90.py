code = r'''
from squircle import *
from methodC import *
from c1 import c1_direct
from mpmath import mp, zeta, mpf, nstr, fabs, matrix, lu_solve
import json
setdps(65)
G = {}
G.update(json.load(open("taugrid34.json")))
G.update(json.load(open("taugrid56.json")))
taus = [0.03, 0.06, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4]
grid = [t for sgn in [1,-1] for t in [sgn*x for x in taus]]
def getZ(tau, s_): return mpf(G[f"{tau},{s_}"])
def lstsq(xs, ys, Q):
    n = len(xs)
    AtA = matrix(Q+1, Q+1); Atb = matrix(Q+1, 1)
    for a in range(Q+1):
        for b in range(Q+1):
            AtA[a,b] = sum(xs[i]**(a+b) for i in range(n))
        Atb[a,0] = sum(xs[i]**a*ys[i] for i in range(n))
    c = lu_solve(AtA, Atb)
    resid = max(abs(ys[i] - sum(c[q,0]*xs[i]**q for q in range(Q+1))) for i in range(n))
    return [c[q,0] for q in range(Q+1)], resid
for s_ in [3,4,5,6]:
    xs = [mpf(t) for t in grid]
    ys = [getZ(t, s_) for t in grid]
    cs, resid = lstsq(xs, ys, 9)
    c1t = c1_direct(s_, K=800)
    print(f"== s={s_} (fit resid {nstr(resid,4)})")
    print("c0:", nstr(cs[0], 22), " vs 4z:", nstr(4*zeta(s_-1), 22), " diff:", nstr(cs[0]-4*zeta(s_-1), 5))
    print("c1:", nstr(cs[1], 22), " vs lat:", nstr(c1t, 22), " diff:", nstr(cs[1]-c1t, 5))
    print("c2:", nstr(cs[2], 22))
    print("c3:", nstr(cs[3], 22))
    print("c4:", nstr(cs[4], 22))
    print("c5:", nstr(cs[5], 22))
'''
exec(code)
