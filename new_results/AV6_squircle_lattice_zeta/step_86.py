code = r'''
from umach2 import eta
from squircle import setdps
from mpmath import mp, nstr, fabs, log, power, mpf, zeta, gamma, matrix, lu_solve, pi
import json
setdps(80)
s = 4
data = {int(k): v for k, v in json.load(open("Uq_data_s4_v2.json")).items()}
def fitU(qi, u0):
    ps = sorted(data.keys())
    xs = [mpf(1)/p_ for p_ in ps]
    ys = [mpf(data[p_][qi])*power(p_, s-1) - u0 for p_ in ps]
    Q = 3
    AtA = matrix(Q, Q); Atb = matrix(Q, 1)
    n = len(xs)
    for a in range(Q):
        for b in range(Q):
            AtA[a,b] = sum(xs[i]**(a+b+2) for i in range(n))
        Atb[a,0] = sum(xs[i]**(a+1)*ys[i] for i in range(n))
    c = lu_solve(AtA, Atb)
    resid = max(abs(ys[i] - sum(c[q,0]*xs[i]**(q+1) for q in range(Q))) for i in range(n))
    return [c[q,0] for q in range(Q)], resid
A3 = zeta(2)*zeta(3) - mpf(29)/16*zeta(5)
u20_pred = zeta(3)*(2*A3)
c20, r20 = fitU(1, u20_pred)
print("u_{2,0} pred:", nstr(u20_pred, 25))
print("u_{2,1} fit:", nstr(c20[0], 20), " u_{2,2} fit:", nstr(c20[1], 18), " u_{2,3}:", nstr(c20[2], 12), " res:", nstr(r20, 4))
# u_{3,0}: fit the constant itself (free)
ps = sorted(data.keys())
xs = [mpf(1)/p_ for p_ in ps]
ys = [mpf(data[p_][2])*power(p_, s-1) for p_ in ps]
Q = 4
AtA = matrix(Q+1, Q+1); Atb = matrix(Q+1, 1)
n = len(xs)
for a in range(Q+1):
    for b in range(Q+1):
        AtA[a,b] = sum(xs[i]**(a+b) for i in range(n))
    Atb[a,0] = sum(xs[i]**a*ys[i] for i in range(n))
c = lu_solve(AtA, Atb)
print("u_{3,0} fit:", nstr(c[0,0], 22), " u_{3,1}:", nstr(c[1,0], 15))
'''
exec(code)
