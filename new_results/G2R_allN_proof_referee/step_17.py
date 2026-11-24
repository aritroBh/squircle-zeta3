import os, json, pickle, shutil, subprocess, time
CWD = "<project>/research/08_theory/G2R_allN_proof_referee"
G2 = "<project>/research/08_theory/G2_M4alpha/"
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
B3 = "<project>/research/10_frontier/W1_corner_family/B_taylor_pi4/"
os.chdir(CWD)
with open("progress.log", "a") as f:
    f.write("step_read_inputs: read result.md, step_05/07/10/11, g2_fit, zsum, hpl, b3_allN_proof, flat_exact, planner_check_b4 -- OK\n")

# set up pristine copy dir
FILES = ["step_10.py", "step_11.py", "g2_fit.py", "mu4_formula.pkl",
         "hpl_fit_k4_(0,0,0).json", "hpl_fit_k4_(1,0,0).json", "b4_pi4_part.json",
         "flat_exact_N640_K4_formula.json"]
for d in ["run_base"]:
    os.makedirs(d, exist_ok=True)
    for fn in FILES:
        shutil.copy2(G2 + fn, os.path.join(d, fn))
print("run_base set up:", os.listdir("run_base"))

# inspect mu4_formula.pkl
terms = pickle.load(open(G2 + "mu4_formula.pkl", "rb"))
print("mu4_formula.pkl:", len(terms), "monomials; sample:", terms[0])
# zeta classification sanity
bad = []
from collections import Counter
zc = Counter()
for mon, cf in terms:
    e2, e3, e4 = mon[8], mon[9], mon[10]
    if not (e2 <= 2 and e3 <= 1 and e4 <= 1): bad.append(mon)
    if sum(1 for x in (e2, e3, e4) if x) > 1: bad.append(mon)   # at most one zeta family
    key = 'z2^2' if e2 == 2 else ('z2' if e2 else ('z3' if e3 else ('z4' if e4 else '1')))
    zc[key] += 1
print("zeta-part monomial counts:", dict(zc), "; inconsistent monomials:", bad)
