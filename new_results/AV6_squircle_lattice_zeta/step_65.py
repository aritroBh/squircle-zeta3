code = r'''
from squircle import setdps
from mpmath import mp, nstr, fabs, log, power, mpf, zeta, gamma, matrix, lu_solve, pi
import json
setdps(100)
s = 4
data = {int(k): v for k, v in json.load(open("Uq_data_s4.json")).items()}
def lstsq(xs, ys, Q, w=None):
    n = len(xs)
    if w is None: w = [1]*n
    AtA = matrix(Q+1, Q+1); Atb = matrix(Q+1, 1)
    for a in range(Q+1):
        for b in range(Q+1):
            AtA[a,b] = sum(w[i]*xs[i]**(a+b) for i in range(n))
        Atb[a,0] = sum(w[i]*xs[i]**a*ys[i] for i in range(n))
    c = lu_solve(AtA, Atb)
    resid = max(abs(ys[i] - sum(c[q,0]*xs[i]**q for q in range(Q+1))) for i in range(n))
    return [c[q,0] for q in range(Q+1)], resid
u1_pred = gamma(s-1)*zeta(s-1)*zeta(s)*(1-power(2,1-s))
A3 = zeta(2)*zeta(3) - mpf(29)/16*zeta(5)
u2_pred = zeta(3)*(2*A3)
ps = sorted(data.keys())
for q in [1,2,3]:
    print(f"== q={q}")
    pred = u1_pred if q==1 else (u2_pred if q==2 else None)
    for psub, Q in [(ps, 4), (ps[2:], 4), (ps[2:], 3), (ps[3:], 3)]:
        xs = [mpf(1)/p_ for p_ in psub]
        ys = [mpf(data[p_][q-1]) * power(p_, s-1) for p_ in psub]
        c, resid = lstsq(xs, ys, Q)
        extra = f" pred-diff={nstr(c[0]-pred, 5)}" if pred else ""
        print(f"  p>={psub[0]} Q={Q}: u0={nstr(c[0], 24)}{extra} (res {nstr(resid,3)})")
        if Q>=3: print(f"        u1={nstr(c[1], 18)}  u2={nstr(c[2], 15)}")
'''
exec(code)
