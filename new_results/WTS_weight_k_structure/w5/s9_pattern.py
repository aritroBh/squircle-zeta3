# Step 9 (supplementary): is k = w-1 (the top weight of M_{w-1}) always the exceptional point of the law?
#   weight 3: law W3(k) = W3(3) - (k-3) S   (k>=3), check k=2;
#   weight 4: quadratic fitted on k=4,5,6, checked on 7,8,9, then compared at k=3.
import sys, json
sys.dont_write_bytecode = True
import mpmath as mp
mp.mp.dps = 140
print("mp.mp.dps =", mp.mp.dps)
ld = lambda w, k: mp.mpf(json.load(open(f"results/w{w}_k{k}_dps140.json"))["value"])
sh = lambda x: mp.nstr(x, 8)
S = mp.zeta(3)/8 + mp.pi**2*mp.log(2)/24 + mp.log(2)**3/12
w3 = {k: ld(3, k) for k in range(2, 12)}
print("weight 3: k>=3 law residuals (k=4..11):", {k: sh(w3[k] - (w3[3] - (k-3)*S)) for k in range(4, 12)})
print("weight 3: k=2 deviation W3(2) - [W3(3) + S] =", mp.nstr(w3[2] - (w3[3] + S), 30))
w4 = {k: ld(4, k) for k in range(3, 10)}
def lag(pk, pv):
    n = len(pk); A = mp.matrix(n, n); y = mp.matrix(n, 1)
    for r, (k, v) in enumerate(zip(pk, pv)):
        for c in range(n): A[r, c] = mp.mpf(k)**c
        y[r] = v
    x = mp.lu_solve(A, y); return [x[i] for i in range(n)]
c = lag([4, 5, 6], [w4[k] for k in (4, 5, 6)])
pv = lambda k: sum(ci*mp.mpf(k)**m for m, ci in enumerate(c))
print("weight 4: quadratic through k=4,5,6; residuals at k=7,8,9:", {k: sh(w4[k]-pv(k)) for k in (7, 8, 9)})
print("weight 4: k=3 deviation W4(3) - P2(3) =", mp.nstr(w4[3]-pv(3), 30))
json.dump(dict(w3_k2_deviation=mp.nstr(w3[2] - (w3[3] + S), 30), w4_k3_deviation=mp.nstr(w4[3]-pv(3), 30),
               w4_quadratic_resid_k7_9={k: sh(w4[k]-pv(k)) for k in (7, 8, 9)}), open("s9_results.json", "w"), indent=1)
