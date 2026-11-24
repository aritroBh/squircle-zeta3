code = r'''
from squircle import setdps
from mpmath import mp, nstr, fabs, log, power, mpf, zeta, gamma, matrix, lu_solve, pi
import json
setdps(80)
s = 4
data = {int(k): v for k, v in json.load(open("Uq_data_s4_v2.json")).items()}
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
u0 = gamma(s-1)*zeta(s-1)*eta(s)
ps = sorted(data.keys())
for psub, Q in [(ps, 4), (ps[1:], 4), (ps[2:], 3)]:
    xs = [mpf(1)/p_ for p_ in psub]
    ys = [mpf(data[p_][0])*power(p_, s-1) - u0 for p_ in psub]
    # fit (D-u0) = u1 x + u2 x^2 + u3 x^3 + u4 x^4 with x = 1/p
    AtA = matrix(Q, Q); Atb = matrix(Q, 1)
    n = len(xs)
    for a in range(Q):
        for b in range(Q):
            AtA[a,b] = sum(xs[i]**(a+b+2) for i in range(n))
        Atb[a,0] = sum(xs[i]**(a+1)*ys[i] for i in range(n))
    c = lu_solve(AtA, Atb)
    resid = max(abs(ys[i] - sum(c[q,0]*xs[i]**(q+1) for q in range(Q))) for i in range(n))
    print(f"pset>={psub[0]}: u1={nstr(c[0,0], 20)}  u2={nstr(c[1,0], 18)}  u3={nstr(c[2,0], 15)}  (res {nstr(resid,3)})")
print("theory u1 = -14.02251936916352983563 (umach)  vs  -13.523...? (part1)")
print("theory u2 = +59.23442463552767339137")
print("theory u3 = -214.7677655889221422249")
'''
exec(code)
