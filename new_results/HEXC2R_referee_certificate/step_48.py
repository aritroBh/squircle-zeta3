
step04 = '''"""step_04.py -- referee TASK 1 (row validity): classify every certificate row by family
(first-emit family in cmzv6.make_relations order: double_shuffle, distribution d=2,3,6,
weight1_branch, lift), and evaluate relation rows NUMERICALLY with MY OWN CMZV evaluator
(reflib truncated multiple sums, 65 dps, cross-validated in step_03): every row used by the
certificates (28+229+1695) plus a random sample of further unused rows."""
import sys, json, time, random
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "../HEXC2_level6_certificate")
sys.path.insert(0, "../HEXC1_hexagon_certificate")
import mpmath as mp
mp.mp.dps = 65
from fractions import Fraction as Fr
import cmzv6
import reflib as R
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

cert23 = json.load(open("../HEXC2_level6_certificate/cert_w23.json"))
cert4 = json.load(open("../HEXC2_level6_certificate/cert_w4.json"))

# --- tag each row with its first-emit family (same emit order as make_relations) ---
def tag_families(n, lower):
    fam = {}
    def emit(row, f):
        r = cmzv6.primitive(row)
        if r and r not in fam: fam[r] = f
    if n == 1:
        for r in cmzv6.weight1_rows():
            rr = cmzv6.primitive(r)
            if rr and rr not in fam: fam[rr] = 'weight1_branch'
    else:
        for a in range(1, n):
            b = n - a
            for u in cmzv6.words(a, False):
                u_adm = u[0] != 0
                if u_adm and a > b: continue
                for v in cmzv6.words(b):
                    if u_adm and a == b and u > v: continue
                    emit(cmzv6.ds(u, v), 'double_shuffle')
    for d in (2, 3, 6):
        for row in cmzv6.distribution(n, d): emit(row, f'distribution_d{d}')
    for k, basis in lower.items():
        for row in basis:
            for w in cmzv6.words(n - k):
                emit(cmzv6.mul(dict(row), {w: 1}), f'lift_w{k}')
    return fam

lower = {}
cert_rows = {2: cert23["2"]["prow"], 3: cert23["3"]["prow"], 4: cert4["prow"]}
random.seed(31415)
summary = {}
for n in (1, 2, 3, 4):
    rows, counts = cmzv6.make_relations(n, lower)
    lower[n] = rows
    if n == 1: continue
    fam = tag_families(n, lower)
    # rows to test: all certificate rows + random unused rows (stratified by family)
    test_idx = list(cert_rows[n])
    unused = [i for i in range(len(rows)) if i not in set(test_idx)]
    random.shuffle(unused)
    per_fam = {}
    for i in unused:
        f = fam[rows[i]]
        if per_fam.get(f, 0) < (12 if n < 4 else 8):
            per_fam[f] = per_fam.get(f, 0) + 1
            test_idx.append(i)
    log(f"weight {n}: testing {len(test_idx)} rows "
        f"({len(cert_rows[n])} cert + {len(test_idx)-len(cert_rows[n])} unused)")
    worst = {}
    ntest = 0
    for i in test_idx:
        row = rows[i]
        for w, c in row:
            assert w[0] != 0 and w[-1] != -1, (n, i, w)
        val = mp.mpc(0)
        for w, c in row:
            val += int(c) * R.cmzv_value(w)
        f = fam[row]
        a = abs(val)
        if a > worst.get(f, (0, 0))[1]: worst[f] = (i, a)
        ntest += 1
    summary[n] = worst
    for f, (i, a) in sorted(worst.items()):
        log(f"    {f:18s} worst |sum| = {mp.nstr(a, 3)}  (row {i})")
    overall = max(a for i, a in worst.values())
    summary[n] = overall
    log(f"    weight {n}: {ntest} rows, overall worst = {mp.nstr(overall, 3)} "
        + ("PASS" if overall < mp.mpf("1e-58") else "FAIL"))
    assert overall < mp.mpf("1e-58")

total = sum(len(cert_rows[n]) for n in (2, 3, 4))
log(f"STEP 04 PASS: all {total} certificate rows + unused samples vanish at 65 dps "
    f"(worst {mp.nstr(max(summary.values()), 3)})")
with open("progress.log", "a") as f:
    f.write("step_04: all %d cert rows + unused samples vanish numerically at 65dps, worst %.1e\\n"
            % (total, max(summary.values())))
'''
with open("steps/step_04.py", "w") as f:
    f.write(step04)
import subprocess
r = subprocess.run(["nice", "-n", "10", "<project>/.venv/bin/python", "steps/step_04.py"], capture_output=True, text=True, timeout=3500)
print(r.stdout)
print(r.stderr[-2500:] if r.returncode else "OK")
