"""verify_g2r.py -- G2R referee verification (< 5 min): the all-n proof run plus one mutant.

Part 1: copy the proof pipeline (step_10.py, step_11.py, g2_fit.py + data) from the read-only
        G2 folder into ./_verify/base and run step_11.py. Expect the final verdict line
        'B4 PROVED for all n (all five parts): True'.
Part 2: MUTATION TEST. In ./_verify/mut, change one rational coefficient of one weight-5 term
        of R4 (c^0 He[0,0,0,0,1]: 1 -> 2) and rerun. The proof MUST fail
        (verdict False, nonempty residual, low-order mismatch).

Exit code 0 iff part 1 passes AND the mutant fails. Runtime ~15 s.
"""
import os, sys, json, re, shutil, subprocess, time

BASE = os.path.dirname(os.path.abspath(__file__))
G2 = "<project>/research/08_theory/G2_M4alpha/"
FILES = ["step_10.py", "step_11.py", "g2_fit.py", "mu4_formula.pkl",
         "hpl_fit_k4_(0,0,0).json", "hpl_fit_k4_(1,0,0).json", "b4_pi4_part.json",
         "flat_exact_N640_K4_formula.json"]

def fresh(dst):
    os.makedirs(dst, exist_ok=True)
    for fn in FILES:
        shutil.copy2(os.path.join(G2, fn), os.path.join(dst, fn))

def run(d):
    t0 = time.time()
    r = subprocess.run([sys.executable, "step_11.py"], cwd=d, capture_output=True, text=True)
    verdict = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else "(no output)"
    return r, verdict, time.time() - t0

t0 = time.time()
ok_all = True

# ---- part 1: pristine ----
d1 = os.path.join(BASE, "_verify", "base")
fresh(d1)
r, verdict, dt = run(d1)
proved = "True" in verdict
ok_all &= proved
print(f"[1] pristine run  ({dt:.1f}s): {verdict}")
for line in r.stdout.splitlines():
    if "even identity" in line or "odd identity" in line or "low orders" in line:
        print("     " + line.strip())
if r.returncode != 0:
    ok_all = False; print("    STDERR:", r.stderr[-400:])

# ---- part 2: mutant (one rational coefficient of one weight-5 term of R4) ----
d2 = os.path.join(BASE, "_verify", "mut")
fresh(d2)
p = os.path.join(d2, "hpl_fit_k4_(0,0,0).json")
data = json.load(open(p))
mutated = False
for i, (nm, v) in enumerate(data):
    m = re.match(r'c\^(-?\d+) He\[(.*)\]', nm)
    if m and len(m.group(2).split(',')) == 5:
        data[i] = [nm, str(int(v) + 1)]
        print(f"[2] mutant: R4 weight-5 term {nm}: coeff {v} -> {int(v)+1}")
        mutated = True
        break
assert mutated
json.dump(data, open(p, "w"))
r, verdict, dt = run(d2)
failed_as_it_should = "False" in verdict
ok_all &= failed_as_it_should
print(f"[2] mutant run  ({dt:.1f}s): {verdict}   (must be False)")
for line in r.stdout.splitlines():
    if "residual" in line or "mismatch" in line:
        print("     " + line.strip())

print(f"VERIFY_G2R: pristine PROVED={proved}, mutant FAILED as required={failed_as_it_should} "
      f"({time.time()-t0:.0f}s total) -> {'PASS' if ok_all else 'FAIL'}")
sys.exit(0 if ok_all else 1)
