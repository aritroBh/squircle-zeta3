import os
work='<project>/research/08_theory/WT4_weight4_all_k'
script = r'''"""WT4 step F: (1) proved closed forms R_D1, R_D2 evaluated at 950 dps vs the 860-digit quadrature values;
(2) k=3 offset: P(3) = W4'' - D1 + D2 vs W_4(3) (quadrature of k=3 kernels AND thm:W4 closed-form eval);
(3) W_4(k) table check k=5..8: formula W4''+(k-4)D1+binom(k-4,2)D2 vs direct quadrature of k-kernels."""
import sys, os, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
import mpmath as mp
from fractions import Fraction as Fr
from qsqrt2 import QS2 as Q, A, B

mp.mp.dps = 950
print("mp.mp.dps =", mp.mp.dps, flush=True)
t0 = time.time()
r2 = mp.sqrt(2); lam = mp.log(2); L = mp.log(1+r2); b = r2-1
pi2 = mp.pi**2; z3 = mp.zeta(3)
Li2b = mp.polylog(2, b); Li3b = mp.polylog(3, b); Li3mb = mp.polylog(3, -b); Li3h = mp.polylog(3, 1/r2)

def Fr_(s): return Fr(s)
lll = json.load(open(os.path.join(HERE, "lll_results.json")))
basis_vals = {
 "log2": lam, "L": L, "pi^2": pi2, "Li2(b)": Li2b, "z3": z3,
 "Li3(b)": Li3b, "Li3(-b)": Li3mb, "Li3(1/sqrt2)": Li3h,
 "Li4(1/2)": mp.polylog(4, mp.mpf('0.5')),
}
def eval_term(nm):
    if nm == "pi^4": return pi2**2
    if nm == "pi^2*log2^2": return pi2*lam**2
    if nm == "log2^4": return lam**4
    parts = nm.split("*")
    v = mp.mpf(1)
    for p in parts:
        v *= basis_vals[p] if not p.startswith("sqrt2") else r2*basis_vals["log2"]  # sqrt2 handled below
    return v
def eval_cf(cf):
    tot = mp.mpf(0)
    for nm, fs in cf.items():
        f = Fr(fs)
        if nm.startswith("sqrt2*"):
            tot += f*r2*eval_term(nm[6:])
        else:
            tot += f*eval_term(nm)
    return tot

R_D1 = eval_cf(lll["D1_closed_form"])
R_D2 = eval_cf(lll["D2_closed_form"])
D1V = mp.mpf(open(os.path.join(HERE, "d1_860.txt")).read().strip())
D2V = mp.mpf(open(os.path.join(HERE, "d2_860.txt")).read().strip())
e1 = abs(R_D1 - D1V); e2 = abs(R_D2 - D2V)
print(f"|R_D1(950dps) - D1(860d)| = {mp.nstr(e1,4)}", flush=True)
print(f"|R_D2(950dps) - D2(860d)| = {mp.nstr(e2,4)}", flush=True)

# ---- W4'' value (M4a 1060-digit quadrature file) ----
W4pp = mp.mpf(open("<project>/research/08_theory/M4a_decomposition/W4pp_1060.txt").read().strip())

# ---- k=3: W_4(3) by quadrature of kernels_deg3[3] at 130 dps ----
kernels = json.load(open(os.path.join(HERE, "kernels_deg3.json")))
def quad_k(k, dps):
    with mp.workdps(dps):
        rr2 = mp.sqrt(2); aN = 1+rr2; bN = rr2-1
        pole = {"0": mp.mpf(0), "a": aN, "-a": -aN, "b": bN, "-b": -bN}
        tot = mp.mpf(0)
        for key, i, j in [("2,1", 2, 1), ("3,0", 3, 0)]:
            ts = [(z, mp.mpf(c[0]) + rr2*mp.mpf(c[1])) for z, c in kernels[str(k)][key]]
            def f(v, ts=ts, i=i, j=j):
                l1 = mp.log((v - bN)*(v + aN)/(4*v))
                l2 = mp.log((aN - v)*(v + bN)/(4*v))
                return l1**i * l2**j * mp.fsum(c/(v - pole[z]) for z, c in ts)
            val = mp.quad(f, [bN, mp.mpf(1)/2, 1, 2, aN], maxdegree=14)
            tot += val
        return +tot
W43_quad = quad_k(3, 130)
print(f"W_4(3) quadrature (130 dps) = {mp.nstr(W43_quad, 60)}", flush=True)

# ---- k=3: W_4(3) via thm:W4 closed form (w4_symbol.conjecture_terms) at 950 dps ----
import w4_symbol as ws
def atom_eval(a):
    if a[0] == 'li':
        z = a[2]
        zm = mp.mpf(str(z.p)) + r2*mp.mpf(str(z.q))
        return mp.polylog(a[1], zm)
    if a[0] == 'log':
        j2, m = ws.logvec(a[1])
        return mp.mpf(float(j2))*lam + mp.mpf(float(m))*L
    if a == ('pi2',): return pi2
    if a == ('z3',): return z3
    raise ValueError(a)
W43_thm = mp.mpf(0)
for c, atoms in ws.conjecture_terms():
    cm = mp.mpf(str(c.p)) + r2*mp.mpf(str(c.q))
    v = cm
    for a in atoms:
        v *= atom_eval(a)
    W43_thm += v
print(f"W_4(3) thm:W4 closed form (950 dps) = {mp.nstr(W43_thm, 60)}", flush=True)
print(f"|W_4(3)quad - W_4(3)thm| = {mp.nstr(abs(W43_quad-W43_thm),4)} (quad at 130 dps)")

# ---- k=3 offset ----
P3 = W4pp - D1V + D2V
delta = P3 - W43_thm
print(f"P(3) = W4'' - D1 + D2 = {mp.nstr(P3, 60)}")
print(f"delta = P(3) - W_4(3) = {mp.nstr(delta, 70)}")
print(f"|delta| > 0: k=3 does NOT lie on the k>=4 quadratic" if abs(delta) > mp.mpf('1e-800') else "k=3 FITS")

# ---- W_4(k) table check k=5..8 ----
from math import comb
tab = {}
for k in range(5, 9):
    pred = W4pp + (k-4)*D1V + comb(k-4, 2)*D2V
    got = quad_k(k, 110)
    e = abs(pred - got)
    tab[k] = {"formula": mp.nstr(pred, 40), "quadrature": mp.nstr(got, 40), "diff": mp.nstr(e, 4)}
    print(f"k={k}: formula={mp.nstr(pred,40)} quad={mp.nstr(got,40)} |diff|={mp.nstr(e,4)} (quad 110 dps)", flush=True)

out = {"R_D1_vs_value": mp.nstr(e1, 6), "R_D2_vs_value": mp.nstr(e2, 6),
       "W43_quad_130": mp.nstr(W43_quad, 130), "W43_thm_950": mp.nstr(W43_thm, 130),
       "W43_quad_minus_thm": mp.nstr(abs(W43_quad-W43_thm), 6),
       "P3": mp.nstr(P3, 130), "delta_P3_minus_W43": mp.nstr(delta, 130),
       "table_check": tab, "dps": mp.mp.dps}
json.dump(out, open(os.path.join(HERE, "sanity_k3_table.json"), "w"), indent=1)
with open(os.path.join(HERE, "progress.log"), "a") as f:
    f.write(f"step_F sanity: R_D1-val={out['R_D1_vs_value']}, R_D2-val={out['R_D2_vs_value']}, delta(k=3)={mp.nstr(delta,8)}, table k=5..8 max diff={max(t['diff'] for t in tab.values())}\n")
print("DONE", time.time()-t0, "s")
'''
p=os.path.join(work,'wt4_F_sanity.py')
open(p,'w').write(script)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','19',PY,'wt4_F_sanity.py'], cwd=work, capture_output=True, text=True, timeout=1700)
print(r.stdout)
print('STDERR:', r.stderr[-1500:])
