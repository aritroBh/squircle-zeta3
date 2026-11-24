code = r'''
from squircle import *
from methodC import *
from mpmath import mp, zeta, mpf, nstr, fabs, log, power, gamma, fac, matrix, lu_solve, pi
import json
setdps(65)
ps_all = [24, 32, 40, 48, 64, 96, 128]
G = {}
G.update(json.load(open("epsgrid34.json")))
G.update(json.load(open("epsgrid56.json")))
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
for s_ in [3,4]:
    u1 = gamma(s_-1)*zeta(s_-1)*zeta(s_)*(1-power(2,1-s_))
    print(f"== s={s_}: predicted d0 = -8s u1 = {nstr(-8*s_*u1, 20)}")
    for pset in [[24,32,40,48,64,96,128],[32,40,48,64,96,128],[48,64,96,128],[64,96,128]]:
        xs = [mpf(1)/p_ for p_ in pset]
        Rs = [ (getZ(p_, s_) - diagpart(mpf(1)/p_, s_))/ (mpf(1)/p_)**s_ for p_ in pset]
        Q = min(3, len(pset)-1)
        d = lstsq(xs, Rs, Q)
        print(f"  pset min={pset[0]}: d0={nstr(d[0], 20)}  d1={nstr(d[1], 12)}  d2={nstr(d[2], 10) if Q>1 else ''}")
'''
exec(code)
