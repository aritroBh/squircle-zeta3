"""w6 step 5: verify the PSLQ closed form of the leading coefficient c_4 and probe the k=5 defect.
All printed numbers come from this script."""
import sys, os, json, glob
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
import mpmath as mp

mp.mp.dps = 110
print("mp.mp.dps =", mp.mp.dps)
lam = mp.log(2); z3 = mp.zeta(3); z5 = mp.zeta(5); pi2 = mp.pi ** 2
c4_closed = (168 * lam ** 6 + 420 * pi2 * lam ** 4 + 378 * pi2 ** 2 * lam ** 2 + 61 * pi2 ** 3
             + 5040 * z3 * lam ** 3 + 2520 * pi2 * z3 * lam + 3780 * z3 ** 2 + 22680 * z5 * lam) / 362880
a110 = json.load(open(os.path.join(W, "analysis_vals_110.json")))
a70 = json.load(open(os.path.join(W, "analysis_vals.json")))
c4_110 = mp.mpf(a110["fit"]["coefficients"][4])
c4_70 = mp.mpf(a70["fit"]["coefficients"][4])
print("c_4 closed form          =", mp.nstr(c4_closed, 100))
print("c_4 fitted (110-dps run) =", mp.nstr(c4_110, 100))
print("closed - fitted(110) =", mp.nstr(c4_closed - c4_110, 5), " (dps", mp.mp.dps, ")")
print("closed - fitted(70)  =", mp.nstr(c4_closed - c4_70, 5), " (dps", mp.mp.dps, ")")

# k=5 defect: delta = W_6(5) - p(5) = -Delta^5 y over k=5..10 (p has degree 4 and fits k>=6 exactly)
y = [mp.mpf(json.load(open(os.path.join(W, "vals_110", f"k{k}_d5.json")))["value"]) for k in range(5, 11)]
d5 = sum((-1) ** (5 - m) * mp.binomial(5, m) * y[m] for m in range(6))
delta = -d5
print("\nk=5 defect delta = W_6(5) - p(5) =", mp.nstr(delta, 60), "(dps", mp.mp.dps, ")")
print("delta / c_4 =", mp.nstr(delta / c4_closed, 30))

mp.mp.dps = 100
Li4 = mp.polylog(4, mp.mpf(1) / 2); Li5 = mp.polylog(5, mp.mpf(1) / 2); Li6 = mp.polylog(6, mp.mpf(1) / 2)
A = [("lam^6", lam ** 6), ("pi^2 lam^4", pi2 * lam ** 4), ("pi^4 lam^2", pi2 ** 2 * lam ** 2), ("pi^6", pi2 ** 3),
     ("z3 lam^3", z3 * lam ** 3), ("pi^2 z3 lam", pi2 * z3 * lam), ("z3^2", z3 ** 2), ("z5 lam", z5 * lam)]
B = A + [("Li4(1/2) lam^2", Li4 * lam ** 2), ("pi^2 Li4(1/2)", pi2 * Li4), ("Li5(1/2) lam", Li5 * lam), ("Li6(1/2)", Li6)]
for name, basis in (("A", A), ("B", B)):
    vec = [delta] + [v for _, v in basis]
    rel = mp.pslq(vec, maxcoeff=10 ** 8, maxsteps=10 ** 6)
    print(f"PSLQ delta vs basis {name} ({len(basis)} elements, dps {mp.mp.dps}):", "no relation (maxcoeff 1e8)" if rel is None else rel)
