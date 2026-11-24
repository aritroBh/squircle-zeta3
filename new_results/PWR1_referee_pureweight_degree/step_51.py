import mpmath as mp, json, sympy as sp
mp.mp.dps = 60
import pwr1_engine as E
w3 = json.load(open('w_exact_n3_tau10.json'))
znum = {E.Z[k]: mp.zeta(k) for k in range(2, 12)}
W3 = {int(k): mp.mpf(str(sp.sympify(v).subs(znum).evalf(50))) for k, v in w3.items()}
vals3 = {mp.mpf(k): mp.mpf(v) for k, v in json.load(open('quad_logS_n3.json')).items()}
def logRat3(tau): return mp.log(1-tau) - mp.log(1-2*tau)
pts = [(t, v - logRat3(t)) for t, v in sorted(vals3.items())]
KMIN, KMAX = 2, 8
A = mp.matrix(len(pts), KMAX-KMIN+1); b = mp.matrix(len(pts), 1)
for i, (t, v) in enumerate(pts):
    for j, k in enumerate(range(KMIN, KMAX+1)): A[i, j] = t**k
    b[i] = v
sol = mp.lu_solve(A.T*A, A.T*b)
print("n=3 pure-part LSQ fit k=2..8 over", len(pts), "points:")
for j, k in enumerate(range(KMIN, KMAX+1)):
    rel = abs(sol[j]-W3[k])/abs(W3[k])
    print(f"  w_{k}(3): fit={mp.nstr(sol[j],25)} relerr={mp.nstr(rel,2)} ({mp.nstr(-mp.log10(rel),3)} digits)")
