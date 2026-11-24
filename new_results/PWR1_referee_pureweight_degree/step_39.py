import mpmath as mp, json
mp.mp.dps = 75
import pwr1_engine as E, sympy as sp

taus = [mp.mpf(x) for x in ['0.002','0.004','0.006','0.008','0.010','0.014','0.020']]
vals = [mp.mpf(s) for s in json.load(open('quad_logP_n2.json'))]

# exact w_k from engine
tot2, pure2, B2 = E.total_log(2)
znum = {E.Z[k]: mp.zeta(k) for k in range(2, 10)}
wex = {}
for k in range(2, 9):
    wex[k] = mp.mpf(sp.N(pure2.coeff(k).subs({s: m for s, m in znum.items()}), 60))
    print(f"exact w_{k} = {mp.nstr(wex[k], 40)}")

# least-squares fit logP(tau) = sum_{k=2}^8 w_k tau^k using subsets; report convergence
def fit(taus, vals, kmax=8):
    A = mp.matrix(len(taus), kmax-1)
    for i, t in enumerate(taus):
        for k in range(2, kmax+1):
            A[i, k-2] = t**k
    b = mp.matrix(vals)
    sol, res = mp.lsq(A, b)  # hmm mpmath lsq? use normal equations via LU
    return sol

# mpmath has no lsq; do interpolation on first 7 points and check residuals on extrapolated points
def interp(taus, vals):
    A = mp.matrix(len(taus), len(taus))
    for i, t in enumerate(taus):
        for j, k in enumerate(range(2, 2+len(taus))):
            A[i, j] = t**k
    b = mp.matrix(len(taus), 1)
    for i, v in enumerate(vals): b[i] = v
    return mp.lu_solve(A, b)

for subset_name, ts, vs in [
    ("first 5 (k=2..6)", taus[:5], vals[:5]),
    ("first 6 (k=2..7)", taus[:6], vals[:6]),
    ("all 7 (k=2..8)", taus, vals),
]:
    sol = interp(ts, vs)
    print(f"\n--- {subset_name}")
    for j, k in enumerate(range(2, 2+len(ts))):
        err = abs(sol[j] - wex[k])
        dig = -mp.log10(err/abs(wex[k])) if err > 0 else mp.inf
        print(f"  w_{k}: fit={mp.nstr(sol[j], 30)}  relerr={mp.nstr(err/abs(wex[k]),2)} ({mp.nstr(dig,3)} digits)")
