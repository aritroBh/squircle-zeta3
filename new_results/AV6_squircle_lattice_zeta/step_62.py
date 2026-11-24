code = r'''
from Uq_direct import *
from squircle import setdps
from mpmath import mp, nstr, fabs, log, power, mpf, zeta, gamma, matrix, lu_solve
import time, json
setdps(80)
s = 4
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
t0=time.time()
ps = [32, 48, 64, 96, 128]
data = {}
for p_ in ps:
    K = 6*p_
    row = []
    for q in [1,2,3]:
        ta=time.time()
        u = Uq_direct(p_, s, q, K=K, R=10)
        row.append(str(u))
        print(f"p={p_} q={q}: {nstr(u, 25)}  ({time.time()-ta:.0f}s)")
    data[p_] = row
json.dump({str(k):v for k,v in data.items()}, open("Uq_data_s4.json","w"))
print("total", time.time()-t0)
'''
exec(code)
with open("progress.log","a") as f: f.write("step_42: U_q(p) direct for s=4, p=32..128\n")
