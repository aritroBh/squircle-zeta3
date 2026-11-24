import mpmath as mp, json, sympy as sp
mp.mp.dps = 100
import pwr1_engine as E

w = json.load(open('w_exact_n2_tau10.json'))
znum = {E.Z[k]: mp.zeta(k) for k in range(2, 12)}
W = {int(k): sp.sympify(v).subs(znum).evalf(60) for k, v in w.items()}
W = {k: mp.mpf(str(v)) for k, v in W.items()}

vals1 = {mp.mpf(k): mp.mpf(v) for k, v in json.load(open('quad_logP_n2_hi.json')).items()}
taus2 = [mp.mpf(x) for x in ['0.002','0.004','0.006','0.008','0.010','0.014','0.020']]
vals2 = dict(zip(taus2, [mp.mpf(s) for s in json.load(open('quad_logP_n2.json'))]))
allv = dict(vals1); allv.update(vals2)
pts = sorted(allv.items())

KMIN, KMAX = 2, 12
A = mp.matrix(len(pts), KMAX-KMIN+1)
b = mp.matrix(len(pts), 1)
for i, (t, v) in enumerate(pts):
    for j, k in enumerate(range(KMIN, KMAX+1)):
        A[i, j] = t**k
    b[i] = v
# least squares via normal equations
sol = mp.lu_solve(A.T*A, A.T*b)
print("least-squares fit k=2..12 over", len(pts), "points; comparison with exact:")
for j, k in enumerate(range(KMIN, KMAX+1)):
    if k in W:
        err = abs(sol[j]-W[k]); rel = err/abs(W[k])
        print(f"  w_{k}: fit={mp.nstr(sol[j],25)} relerr={mp.nstr(rel,2)} ({mp.nstr(-mp.log10(rel),3)} digits)")
    else:
        print(f"  w_{k}: fit={mp.nstr(sol[j],25)} (no exact ref; D16 json: see below)")
json.dump([str(sol[j]) for j in range(KMAX-KMIN+1)], open('fit_w2_w12.json','w'))
