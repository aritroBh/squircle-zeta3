code = r'''
from Uq_direct import *
from squircle import setdps
from mpmath import mp, nstr, fabs, log, power, mpf, zeta, gamma, matrix, lu_solve, pi
import json
setdps(80)
s = 4
data = {int(k): v for k, v in json.load(open("Uq_data_s4.json")).items()}
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
ps = sorted(data.keys())
xs = [mpf(1)/p_ for p_ in ps]
print("s = 4")
# u_1 predictions
u1_pred = gamma(s-1)*zeta(s-1)*zeta(s)*(1-power(2,1-s))
A3 = zeta(2)*zeta(3) - mpf(29)/16*zeta(5)
u2_pred = zeta(3)*(2*A3)  # zeta(s-1) L_2(s-1), L_2(3) = 2 A(3)
print("u1 pred:", nstr(u1_pred, 25))
print("u2 pred:", nstr(u2_pred, 25), " A(3) =", nstr(A3, 25))
fits = {}
for qi, q in enumerate([1,2,3]):
    ys = [mpf(data[p_][qi]) * power(p_, s-1) for p_ in ps]
    c = lstsq(xs, ys, 4)
    fits[q] = c
    print(f"q={q}: u0={nstr(c[0], 25)}")
    print(f"      u1={nstr(c[1], 22)}")
    print(f"      u2={nstr(c[2], 18)}")
    print(f"      u3={nstr(c[3], 14)}")
print("u_1^(0) diff:", nstr(fits[1][0]-u1_pred, 6))
print("u_2^(0) diff:", nstr(fits[2][0]-u2_pred, 6))
'''
exec(code)
