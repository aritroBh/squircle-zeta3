"""verify_files.py  -- cheap consistency checks on the author's STORED numbers (read only).
 1. values quoted in the compute-step report vs results/w4_k<k>_dps120.json
 2. third differences (k=3..8, dps 90) and the quadratic fit through k=4,5,6 (dps 120) recomputed from the stored values;
    coefficients compared with the ones quoted in the report
 3. accuracy of the dps-90 / dps-120 values: compare with the dps-360 runs (k=3..8) and with the dps-360 fit (k=9,10)
PRECISION: mp.mp.dps set after all imports and printed with every residual."""
import sys, os, json
sys.dont_write_bytecode = True
import mpmath as mp

W4 = "<project>/research/08_theory/WTS_weight_k_structure/w4"
R = os.path.join(W4, "results")
mp.mp.dps = 400                                   # after all imports
print("mp.mp.dps =", mp.mp.dps)


def val(k, d):
    return mp.mpf(json.load(open(os.path.join(R, f"w4_k{k}_dps{d}.json")))["value"])


QUOTED = {3: "0.79578505983224081132035981962972676254652741387102124725072771573639123156493438773599028533900548032786",
          10: "-17.240722110716741024925885849398081113156490982529840283182275241420768455633637530814078397642729187265",
          6: "-0.7000559105407769521880758201067125626994755610772090774082125076457630774606425701474515751126549729753"}
for k, s in QUOTED.items():
    print(f"k={k}: quoted value - stored dps120 value = {mp.nstr(mp.mpf(s) - val(k, 120), 5)}  (quoted string has ~{len(s)-3} digits; mp.dps {mp.mp.dps})")

v90 = {k: val(k, 90) for k in range(3, 11)}
v120 = {k: val(k, 120) for k in range(3, 11)}
d3 = [v90[k + 3] - 3 * v90[k + 2] + 3 * v90[k + 1] - v90[k] for k in range(3, 8)]
print("third differences (dps-90 stored) starting at k=3..7:", [mp.nstr(x, 6) for x in d3])
d3b = [v120[k + 3] - 3 * v120[k + 2] + 3 * v120[k + 1] - v120[k] for k in range(3, 8)]
print("third differences (dps-120 stored) starting at k=3..7:", [mp.nstr(x, 6) for x in d3b])


def fit(vals, kk):
    A = mp.matrix(3, 3); y = mp.matrix(3, 1)
    for r, k in enumerate(kk):
        for c in range(3):
            A[r, c] = mp.mpf(k) ** c
        y[r] = vals[k]
    return mp.lu_solve(A, y)


c = fit(v120, [4, 5, 6])
QC = ["-6.6181741779684916564133942925781656214736139387661725375138737447690932225663943745728226570527894904735",
      "4.05926480134045186537075643038609276490832254726828493689428573996580551208748305861342631189688033680965",
      "-0.5121519594615276802222005586068084314076610251644651711461125889630973035394207374237551885955874306488"]
for j in range(3):
    print(f"c_{j}: my fit through k=4,5,6 of the stored dps-120 values - quoted = {mp.nstr(c[j] - mp.mpf(QC[j]), 5)}  (mp.dps {mp.mp.dps})")
print("fit(4,5,6) reproduces k=7..10 (stored dps120):",
      [mp.nstr(sum(c[j] * mp.mpf(k) ** j for j in range(3)) - v120[k], 4) for k in range(7, 11)])

# accuracy of stored values against dps-360 runs (k=3..8)
for d in (90, 120):
    worst = max(abs(val(k, d) - val(k, 360)) for k in range(3, 9))
    print(f"max_k=3..8 |stored dps{d} value - stored dps360 value| = {mp.nstr(worst, 5)}  (so ~{int(-mp.log10(worst))} digits right)")
c360 = fit({k: val(k, 360) for k in range(3, 9)}, [4, 5, 6])
for k in (9, 10):
    p = sum(c360[j] * mp.mpf(k) ** j for j in range(3))
    print(f"k={k}: stored dps120 value - quadratic from the dps-360 fit = {mp.nstr(val(k, 120) - p, 5)}   dps90: {mp.nstr(val(k, 90) - p, 5)}")
print("stored-vs-quoted c_j: c_j(dps360 fit) - c_j(dps120 fit) =", [mp.nstr(c360[j] - c[j], 4) for j in range(3)])
