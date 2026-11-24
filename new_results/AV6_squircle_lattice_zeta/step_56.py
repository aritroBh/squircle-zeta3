code = r'''
from squircle import *
from methodC import *
from mpmath import mp, zeta, mpf, nstr, fabs, log, power, gamma, fac, matrix, lu_solve, pi
import json
setdps(65)
G = json.load(open("epsgrid_all.json"))
def getZ(p_, s_): return mpf(G[f"{p_},{s_}"])
def diagpart(eps, s_):
    return 8*zeta(s_-1) - 4*zeta(s_)*(1-power(2,-s_*eps))
def lstsq(xs, ys, Q):
    n = len(xs)
    A = matrix(n, Q+1)
    for i in range(n):
        xp = mpf(1)
        for q in range(Q+1):
            A[i,q] = xp; xp *= xs[i]
    AtA = matrix(Q+1, Q+1); Atb = matrix(Q+1, 1)
    for a in range(Q+1):
        for b in range(Q+1):
            AtA[a,b] = sum(A[i,a]*A[i,b] for i in range(n))
        Atb[a,0] = sum(A[i,a]*ys[i] for i in range(n))
    c = lu_solve(AtA, Atb)
    return [c[q,0] for q in range(Q+1)]
for s_ in [3,4,5,6]:
    print(f"== s={s_}")
    for pset, Q in [([48,64,96,128,160,192,256],5),([64,96,128,160,192,256],4),([96,128,160,192,256],3)]:
        xs = [mpf(1)/p_ for p_ in pset]
        Rs = [(getZ(p_, s_) - diagpart(mpf(1)/p_, s_))/(mpf(1)/p_)**s_ for p_ in pset]
        d = lstsq(xs, Rs, Q)
        print(f"  >={pset[0]} Q={Q}: d1={nstr(d[1], 16)}  d2={nstr(d[2], 14)}  d3={nstr(d[3], 12) if Q>2 else ''}")
'''
exec(code)
