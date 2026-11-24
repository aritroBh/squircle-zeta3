# Step 5: polynomial-law analysis for W_5(k) (values from s2_eval.py at dps 140 and 110).
#   (a) accuracy of the stored values (dps110 vs dps140; W5(4) vs the 600-digit reference),
#   (b) the REQUESTED test: smallest d with vanishing (d+1)-th finite difference on k=4..9, fit through first d+1 points,
#   (c) the cubic on k>=5: fit through k=5..8, check k=9,10,11; direct c_m by integrating the exact cubic kernels,
#   (d) the k=4 defect, (e) PSLQ for the leading coefficient and the defect.
import sys, os, json, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import fastred                      # (only for unpickling Fractions-based data; no precision side effects)
import mpmath as mp
import wquad
mp.mp.dps = 140                     # set AFTER all imports
print("mp.mp.dps =", mp.mp.dps)
def load(k, dps): return mp.mpf(json.load(open(f"results/w5_k{k}_dps{dps}.json"))["value"])
def sh(x, n=6): return mp.nstr(x, n)
ks = list(range(4, 12))
V = {k: load(k, 140) for k in ks}; V110 = {k: load(k, 110) for k in ks}
res = {}
print("(a) |W5(k;dps140) - W5(k;dps110)| =", {k: sh(V[k]-V110[k], 3) for k in ks})
ref = mp.mpf(open("../../M4a_decomposition/W5_600.txt").read().strip())
print("    W5(4;dps140) - W5_600.txt =", sh(V[4]-ref, 4), " [mp.mp.dps =", mp.mp.dps, "]")

def fd(vals, m):
    s = list(vals)
    for _ in range(m): s = [s[i+1]-s[i] for i in range(len(s)-1)]
    return s
def lagfit(pk, pv):
    n = len(pk)
    A = mp.matrix(n, n); y = mp.matrix(n, 1)
    for r, (k, v) in enumerate(zip(pk, pv)):
        for c in range(n): A[r, c] = mp.mpf(k)**c
        y[r] = v
    x = mp.lu_solve(A, y)
    return [x[i] for i in range(n)]
def pv(cs, k): return sum(c*mp.mpf(k)**m for m, c in enumerate(cs))

print("(b) REQUESTED window k=4..9:")
w = list(range(4, 10))
for m in range(1, 6):
    print(f"    Delta^{m} = " + ", ".join(sh(x, 8) for x in fd([V[k] for k in w], m)))
for dd in (0, 1, 2, 3, 4):
    d1 = fd([V[k] for k in w], dd+1)
    print(f"    d={dd}: (d+1)-th differences max|.| = {sh(max(abs(x) for x in d1), 6)}")
cs4 = lagfit([4, 5, 6, 7], [V[k] for k in (4, 5, 6, 7)])
print("    cubic through k=4,5,6,7: residual at k=8:", sh(V[8]-pv(cs4, 8), 8), " at k=9:", sh(V[9]-pv(cs4, 9), 8))
cs4q = lagfit([4, 5, 6, 7, 8], [V[k] for k in (4, 5, 6, 7, 8)])
print("    quartic through k=4..8: residual at k=9:", sh(V[9]-pv(cs4q, 9), 8), " at k=10:", sh(V[10]-pv(cs4q, 10), 8))
res["requested_window"] = dict(delta4=[sh(x, 12) for x in fd([V[k] for k in w], 4)], delta5=[sh(x, 12) for x in fd([V[k] for k in w], 5)],
                               cubic_through_4_7_resid_k8=sh(V[8]-pv(cs4, 8), 10), cubic_through_4_7_resid_k9=sh(V[9]-pv(cs4, 9), 10))

print("(c) cubic on k>=5:")
cs5 = lagfit([5, 6, 7, 8], [V[k] for k in (5, 6, 7, 8)])
for k in (9, 10, 11):
    print(f"    fit through k=5..8: residual at k={k}: {sh(V[k]-pv(cs5, k), 4)}")
for m in range(4, 7):
    print(f"    Delta^{m} on k=5..11: " + ", ".join(sh(x, 4) for x in fd([V[k] for k in range(5, 12)], m)))
res["cubic_k_ge_5"] = dict(resid_k9=sh(V[9]-pv(cs5, 9), 4), resid_k10=sh(V[10]-pv(cs5, 10), 4), resid_k11=sh(V[11]-pv(cs5, 11), 4),
                           delta4_k5_11=[sh(x, 4) for x in fd([V[k] for k in range(5, 12)], 4)])

# direct kernel integration of the exact cubic kernel
K = pickle.load(open("exact_cubic_kernel_w5.pkl", "rb"))
sq2 = mp.sqrt(2)
def kernel_terms(sel):    # sel(m) -> function giving the (key,pole,comp) coefficient
    acc = {}
    for (key, pole, comp), cs in K["cub"].items():
        c = sel(cs)
        acc.setdefault((key, pole), mp.mpf(0))
        val = mp.mpf(c.numerator)/c.denominator
        acc[(key, pole)] += val if comp == "p0" else val*sq2
    return [(key[0], key[1], pole, c) for (key, pole), c in acc.items() if c != 0]
cm = []
for m in range(4):
    terms = kernel_terms(lambda cs, m=m: cs[m])
    val, err = wquad.integrate_terms(terms)
    cm.append(val)
    print(f"    c_{m} (direct kernel integral) = {mp.nstr(val, 50)}   |c_{m} - fit| = {sh(val-cs5[m], 4)}")
res["cm_direct"] = [mp.nstr(c, 135) for c in cm]
res["cm_fit"] = [mp.nstr(c, 135) for c in cs5]
res["cm_fit_minus_direct"] = [sh(cs5[m]-cm[m], 4) for m in range(4)]

# defect at k=4
dterms = {}
for (key, pole, comp), dv in K["defect"].items():
    val = mp.mpf(dv.numerator)/dv.denominator
    dterms[(key, pole)] = dterms.get((key, pole), mp.mpf(0)) + (val if comp == "p0" else val*sq2)
dval, _ = wquad.integrate_terms([(key[0], key[1], pole, c) for (key, pole), c in dterms.items()])
delta_fd = V[4] - pv(cs5, 4)
print("(d) k=4 defect: W5(4) - P_cubic(4) =", mp.nstr(delta_fd, 50))
print("    direct integral of defect kernel =", mp.nstr(dval, 50), " diff:", sh(delta_fd-dval, 4))
res["defect_k4"] = mp.nstr(delta_fd, 135); res["defect_k4_direct_minus_fd"] = sh(delta_fd-dval, 4)

# PSLQ
lam = mp.log(2); L = mp.log(1+sq2); pi = mp.pi; z3 = mp.zeta(3); z5 = mp.zeta(5)
def pslq(x, names, vals, tol=mp.mpf(10)**-100, maxc=10**8):
    r = mp.pslq([x] + vals, tol=tol, maxcoeff=maxc, maxsteps=10**6)
    print(f"    PSLQ({mp.nstr(x,20)}) vs {names}: {r}")
    return r
print("(e) PSLQ")
bas_names = ["z5", "pi^2 z3", "pi^4 lam", "z3 lam^2", "pi^2 lam^3", "lam^5"]
bas = [z5, pi**2*z3, pi**4*lam, z3*lam**2, pi**2*lam**3, lam**5]
res["pslq_c3_lam_basis"] = pslq(cm[3], bas_names, bas)
bas2_names = bas_names + ["pi^4 L", "z3 L^2", "pi^2 L^3", "L^5", "z3 lam L", "pi^2 lam L^2", "lam^2 L^3", "lam^3 L^2","lam^4 L","lam L^4","pi^2 lam^2 L","pi^4 L"]
json.dump(res, open("s5_results.json", "w"), indent=1, default=str)
