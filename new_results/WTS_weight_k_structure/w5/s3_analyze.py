# Step 3: controls (weight 2, weight 3, W5(4) reference) and finite-difference / polynomial analysis.
import sys, os, json, glob
sys.dont_write_bytecode = True
import mpmath as mp
DPS = 140
mp.mp.dps = DPS            # no module here resets it; printed below next to every residual
def load(w, k, dps):
    return mp.mpf(json.load(open(f"results/w{w}_k{k}_dps{dps}.json"))["value"])
def sh(x): return mp.nstr(x, 6)
sq2 = mp.sqrt(2); lam = mp.log(2); L = mp.log(1 + sq2); b = sq2 - 1; pi = mp.pi; z3 = mp.zeta(3)
print("mp.mp.dps =", mp.mp.dps)
out = {"controls": [], "dps": DPS}
# ---- weight 2 ----
W2 = (pi**2 + 6*lam**2) / 24
for k in (2, 3, 5, 8):
    r = load(2, k, 140) - W2   # stored at dps 80
    print(f"W2 control k={k}: W_2(k) - (pi^2+6 log^2 2)/24 = {sh(r)}   [mp.mp.dps={mp.mp.dps}, data computed at 140 digits]")
    out["controls"].append(dict(name=f"weight2_k{k}", residual=sh(r), pass_=bool(abs(r) < mp.mpf(10)**-130)))
# ---- weight 3 ----
Li3 = lambda x: mp.polylog(3, x); Li2 = lambda x: mp.polylog(2, x)
W3p = ((7*sq2/32 - mp.mpf(3)/16)*z3 + (1+sq2)*Li3(b) + (1-sq2)*Li3(-b) - 2*sq2*Li3(1/sq2)
       + (2*L - sq2*lam)*Li2(b) + (sq2/12 - mp.mpf(1)/8)*pi**2*L - (1+sq2)/24*pi**2*lam
       + (1+sq2)/24*lam**3 + sq2/4*lam**2*L - (1+sq2)/2*lam*L**2 + (5-sq2)/6*L**3)
S = z3/8 + pi**2*lam/24 + lam**3/12
for k in range(3, 12):
    r = load(3, k, 140) - (W3p - (k-3)*S)
    print(f"W3 control k={k}: W_3(k) - [W3' - (k-3) S] = {sh(r)}   [mp.mp.dps={mp.mp.dps}, data at 140 digits]")
    out["controls"].append(dict(name=f"weight3_k{k}", residual=sh(r), pass_=bool(abs(r) < mp.mpf(10)**-130)))
# ---- W5(4) reference ----
ref = mp.mpf(open("../../M4a_decomposition/W5_600.txt").read().strip())
for dps in (80, 110, 140):
    r = load(5, 4, dps) - ref
    print(f"W5(4) control (data dps {dps}): W_5(4) - W5_600.txt = {sh(r)}   [mp.mp.dps={mp.mp.dps}]")
    out["controls"].append(dict(name=f"W5_k4_vs_W5_600_dps{dps}", residual=sh(r), pass_=bool(abs(r) < mp.mpf(10)**-(dps-5))))
# ---- cross precision ----
ks = list(range(4, 12))
v80 = {k: load(5, k, 80) for k in ks}; v110 = {k: load(5, k, 140) for k in ks}   # v110 name kept: holds dps-140 data here
print("cross-precision |W5(k;dps80) - W5(k;dps140)|:", {k: sh(v80[k]-v110[k]) for k in ks})
# ---- finite differences ----
def diffs(vals):
    rows = [list(vals)]
    while len(rows[-1]) > 1:
        r = rows[-1]; rows.append([r[i+1]-r[i] for i in range(len(r)-1)])
    return rows
def report(window, tag):
    rows = diffs([v110[k] for k in window])
    print(f"--- finite differences, window k={window[0]}..{window[-1]} ({tag}), mp.mp.dps={mp.mp.dps}")
    for m, r in enumerate(rows):
        print(f"  Delta^{m}: " + ", ".join(mp.nstr(x, 12) for x in r))
    return rows
rowsA = report(list(range(4, 10)), "REQUESTED window")
rowsB = report(list(range(4, 12)), "extended, includes k=4")
rowsC = report(list(range(5, 12)), "extended, EXCLUDES k=4")
out["cross_precision_80_vs_140"] = {k: sh(v80[k]-v110[k]) for k in ks}
json.dump(out, open("s3_controls.json", "w"), indent=1)
