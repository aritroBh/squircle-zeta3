import importlib, pwr1_engine as E
importlib.reload(E)
import sympy as sp, json

# n=2: total log and pure part
tot, pure, B = E.total_log(2)
print("n=2: log(S/S(1)) = logRat + pure; Rat_2 should be 1 (logRat=0):")
print("logRat_2 =", E.logRat(2).d)
print("\npure part coefficients w_k(2) = [tau^k] log(P/(4 sqrt2)):")
for k in range(1, E.NMAX+1):
    c = pure.coeff(k)
    gam = sp.expand(c).coeff(E.GAM, 1)
    ok_w = E.weight_check(c, k)
    print(f"k={k}: {c}   [GAM-part={gam}, weight-{k} homogeneous: {ok_w}]")
