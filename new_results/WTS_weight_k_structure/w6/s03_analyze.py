"""w6 step 3: controls, finite differences, exact polynomial fit.

usage: s03_analyze.py [VALSDIR] [DPS]      (default vals, 70)
All printed numbers come from this script.  Writes analysis_<VALSDIR>.json.
"""
import sys, os, json, glob
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
from fractions import Fraction as Fr
import mpmath as mp

VALS = sys.argv[1] if len(sys.argv) > 1 else "vals"
DPS = int(sys.argv[2]) if len(sys.argv) > 2 else 70
mp.mp.dps = DPS          # AFTER all imports
print(f"VALS={VALS}  mp.mp.dps = {mp.mp.dps}")

# ---------------------------------------------------------------- load data
vals = {}      # (k,d) -> mpf
strs = {}
errs = {}
for fn in glob.glob(os.path.join(W, VALS, "k*_d*.json")):
    r = json.load(open(fn))
    vals[(r["k"], r["d"])] = mp.mpf(r["value"])
    strs[(r["k"], r["d"])] = r["value"]
    errs[(r["k"], r["d"])] = mp.mpf(r["err_est"])
    assert r["dps"] == DPS, (fn, r["dps"])
ks6 = sorted(k for (k, d) in vals if d == 5)
print("weight-6 data at k =", ks6)
print("max quadrature err_est over all stored pieces:", mp.nstr(max(errs.values()), 3))

# ---------------------------------------------------------------- constants
lam = mp.log(2)
r2 = mp.sqrt(2)
b = r2 - 1
L = mp.log(1 + r2)
z3 = mp.zeta(3)
pi2 = mp.pi ** 2
W2exact = (pi2 + 6 * lam ** 2) / 24
S3 = z3 / 8 + pi2 * lam / 24 + lam ** 3 / 12
# thm:M3 of paper/squircle_note/main.tex, transcribed from the LaTeX
W3p = ((7 * r2 / 32 - mp.mpf(3) / 16) * z3 + (1 + r2) * mp.polylog(3, b) + (1 - r2) * mp.polylog(3, -b)
       - 2 * r2 * mp.polylog(3, 1 / r2) + (2 * L - r2 * lam) * mp.polylog(2, b)
       + (r2 / 12 - mp.mpf(1) / 8) * pi2 * L - (1 + r2) / 24 * pi2 * lam + (1 + r2) / 24 * lam ** 3
       + r2 / 4 * lam ** 2 * L - (1 + r2) / 2 * lam * L ** 2 + (5 - r2) / 6 * L ** 3)

controls = []
def ctl(name, resid, tol=mp.mpf(10) ** (-(DPS - 12))):
    ok = abs(resid) < tol
    controls.append({"name": name, "residual": mp.nstr(resid, 5), "pass": bool(ok), "dps": mp.mp.dps})
    print(f"CONTROL {'PASS' if ok else 'FAIL'}  {name}: residual {mp.nstr(resid, 5)} (dps {mp.mp.dps}, tol {mp.nstr(tol, 2)})")

# ---- control A: weight 2, all k
for k in sorted(k for (k, d) in vals if d == 1):
    ctl(f"W_2(k={k}) - (pi^2+6 log^2 2)/24", vals[(k, 1)] - W2exact)

# ---- control B: weight 3 = W3' - (k-3) S, k >= 3
for k in sorted(k for (k, d) in vals if d == 2):
    ctl(f"W_3(k={k}) - [W3' - (k-3)*S]", vals[(k, 2)] - (W3p - (k - 3) * S3))

# ---- control C: k=3 sanity against the paper's printed digits (M3 = W4 + W3' + W2)
if (3, 3) in vals:
    W4_paper = mp.mpf("0.795785059832240811320359819629726762546527413871")
    M3_paper = mp.mpf("0.022770521491786577456730509742327800250183691729217")
    ctl("W_4(k=3) - W_4 printed in thm:M3 (48 digits)", vals[(3, 3)] - W4_paper, tol=mp.mpf(10) ** -46)
    M3_mine = vals[(3, 3)] + vals[(3, 2)] + vals[(3, 1)]
    ctl("W_4+W_3+W_2 at k=3 - M3 printed in thm:M3 (51 digits)", M3_mine - M3_paper, tol=mp.mpf(10) ** -49)

# ---- control D: k=5 against the stored independent slow-reducer / direct-integrand run (M5a folder)
M5A = os.path.join(W, "..", "..", "M5a_next_coefficient")
q5 = json.load(open(os.path.join(M5A, "quad5_100.json")))
stored = {5: mp.mpf(q5["W6"]), 4: mp.mpf(q5["W5pp"]), 3: mp.mpf(q5["W4ppp"]), 2: mp.mpf(q5["W35_quad"]), 1: mp.mpf(q5["T"])}
for d in (5, 4, 3, 2, 1):
    if (5, d) in vals:
        ctl(f"W_{d+1}(k=5) - stored M5a quad5_100.json (100 digits)", vals[(5, d)] - stored[d], tol=mp.mpf(10) ** -(min(DPS, 100) - 12))
if all((5, d) in vals for d in (1, 2, 3, 4, 5)):
    M5mine = sum(vals[(5, d)] for d in (1, 2, 3, 4, 5))
    M5txt = mp.mpf(open(os.path.join(M5A, "M5_400.txt")).read().strip())
    ctl("sum_{w=2..6} W_w(5) - M5 (M5a/M5_400.txt, direct integrand quadrature)", M5mine - M5txt,
        tol=mp.mpf(10) ** -(min(DPS, 100) - 12))
    print("   M5 (stored) =", mp.nstr(M5txt, 40), "  sum of my parts =", mp.nstr(M5mine, 40))

# ---------------------------------------------------------------- finite differences
def diffs(seq, n):
    cur = list(seq)
    for _ in range(n):
        cur = [cur[i + 1] - cur[i] for i in range(len(cur) - 1)]
    return cur

def fd_report(label, ks, ys, nmax=9):
    print(f"--- finite differences of {label}, k = {ks[0]}..{ks[-1]} ({len(ks)} points)")
    res = {}
    for n in range(1, min(nmax, len(ys) - 1) + 1):
        dn = diffs(ys, n)
        m = max(abs(x) for x in dn)
        res[n] = (m, dn)
        print(f"   order {n}: {len(dn)} windows, max |Delta^{n}| = {mp.nstr(m, 8)}; first = {mp.nstr(dn[0], 25)}")
    return res

thresh = mp.mpf(10) ** (-(DPS - 18))
def smallest_degree(ks, ys):
    for n in range(1, len(ys)):
        dn = diffs(ys, n)
        if max(abs(x) for x in dn) < thresh:
            return n - 1, max(abs(x) for x in dn)
    return -1, None

out = {"vals_dir": VALS, "dps": DPS, "threshold": mp.nstr(thresh, 3), "controls": controls}
req = [k for k in range(5, 11) if k in ks6]
sets = {"requested_k5_10": req}
if all(k in ks6 for k in range(6, 11)):
    sets["k6_10"] = list(range(6, 11))
ext_all = [k for k in range(5, 16) if k in ks6]
if len(ext_all) > len(req):
    sets["all_k5_to_%d" % ext_all[-1]] = ext_all
    sets["k6_to_%d" % ext_all[-1]] = [k for k in ext_all if k >= 6]
out["degree_search"] = {}
for name, kk in sets.items():
    ys = [vals[(k, 5)] for k in kk]
    print(f"\n=== set {name}: k = {kk}")
    fd_report("W_6(k)", kk, ys)
    d, m = smallest_degree(kk, ys)
    print(f"   smallest d with (d+1)-th difference < {mp.nstr(thresh,2)}: d = {d}" + (f" (residual {mp.nstr(m,5)})" if m is not None else " (none found)"))
    out["degree_search"][name] = {"k": kk, "degree": d, "residual": (mp.nstr(m, 5) if m is not None else None),
        "max_abs_diff_by_order": {str(n): mp.nstr(max(abs(x) for x in diffs(ys, n)), 6) for n in range(1, len(ys))}}

# ---------------------------------------------------------------- exact Lagrange fit (rational arithmetic)
def to_fr(k):
    return Fr(strs[(k, 5)])     # exact decimal expansion of the stored value string

def poly_through(pts):
    """monomial coefficients (low to high, Fractions) of the interpolating polynomial through (k, y) Fractions"""
    n = len(pts)
    coeffs = [Fr(0)] * n
    for a in range(n):
        ka, ya = pts[a]
        num = [Fr(1)]
        den = Fr(1)
        for c in range(n):
            if c == a: continue
            kc = pts[c][0]
            new = [Fr(0)] * (len(num) + 1)
            for t, v in enumerate(num):
                new[t + 1] += v
                new[t] -= kc * v
            num = new
            den *= (ka - kc)
        for t, v in enumerate(num):
            coeffs[t] += ya * v / den
    return coeffs

def fr_to_mp(f):
    return mp.mpf(f.numerator) / mp.mpf(f.denominator)

def fit(kfit, label):
    pts = [(Fr(k), to_fr(k)) for k in kfit]
    c = poly_through(pts)
    return c

d_final = None
fits = {}
for name in ("k6_to_%d" % ext_all[-1], "k6_10", "requested_k5_10"):
    if name in out["degree_search"] and out["degree_search"][name]["degree"] >= 0:
        d_final = out["degree_search"][name]["degree"]
        regime = name
        break
if d_final is not None:
    kk = out["degree_search"][regime]["k"]
    print(f"\n=== exact fit, degree d = {d_final}, regime used: {regime} (k = {kk})")
    kfit = kk[: d_final + 1]
    c = fit(kfit, "first")
    print("fit through k =", kfit)
    for j, cj in enumerate(c):
        print(f"   c_{j} = {mp.nstr(fr_to_mp(cj), DPS - 5)}")
    # remaining points
    fitres = {}
    for k in sorted(k for k in ks6):
        pk = sum(c[j] * Fr(k) ** j for j in range(len(c)))
        r = vals[(k, 5)] - fr_to_mp(pk)
        fitres[k] = r
        tag = "fitted" if k in kfit else ("check" if k in kk else "outside regime")
        print(f"   k={k:2d} [{tag}]: W_6(k) - poly(k) = {mp.nstr(r, 6)}   (dps {mp.mp.dps})")
    # alternative fit from the LAST d+1 points of the regime, for coefficient stability
    kfit2 = kk[-(d_final + 1):]
    c2 = fit(kfit2, "last")
    print("alternative fit through k =", kfit2)
    relmax = mp.mpf(0)
    for j in range(len(c)):
        diff = abs(fr_to_mp(c[j]) - fr_to_mp(c2[j]))
        rel = diff / max(mp.mpf(1), abs(fr_to_mp(c[j])))
        relmax = max(relmax, rel)
        print(f"   c_{j}: |first-fit - last-fit| = {mp.nstr(diff, 5)}  (relative to max(1,|c|): {mp.nstr(rel, 5)})")
    print("   worst relative difference between the two fits:", mp.nstr(relmax, 5))
    out["fit"] = {"regime": regime, "regime_k": kk, "degree": d_final, "fit_points": kfit,
                  "coefficients": [mp.nstr(fr_to_mp(cj), DPS - 5) for cj in c],
                  "coefficients_last_fit": [mp.nstr(fr_to_mp(cj), DPS - 5) for cj in c2],
                  "worst_relative_difference_first_vs_last_fit": mp.nstr(relmax, 5),
                  "residuals": {str(k): mp.nstr(r, 6) for k, r in fitres.items()},
                  "k5_prediction": mp.nstr(fr_to_mp(sum(c[j] * Fr(5) ** j for j in range(len(c)))), 30),
                  "k5_actual": mp.nstr(vals[(5, 5)], 30)}
    print("   k=5 polynomial prediction:", out["fit"]["k5_prediction"], " actual:", out["fit"]["k5_actual"])
    # lower weights, same procedure on the data we have (bonus structure check)
    for d in (4, 3, 2):
        kk2 = sorted(k for (k, dd) in vals if dd == d)
        ys2 = [vals[(k, d)] for k in kk2]
        reg = [k for k in kk2 if k >= d + 1]
        ys3 = [vals[(k, d)] for k in reg]
        dg, m = smallest_degree(reg, ys3)
        print(f"   [bonus] weight {d+1}: k = {kk2[0]}..{kk2[-1]}; regime k>={d+1}: smallest degree = {dg}" + (f" (resid {mp.nstr(m,4)})" if m is not None else " (none or too few points)"))
        out.setdefault("bonus_lower_weights", {})[str(d + 1)] = {"regime_k": reg, "degree": dg}
json.dump(out, open(os.path.join(W, f"analysis_{VALS}.json"), "w"), indent=1)
print("\nsaved", f"analysis_{VALS}.json")
