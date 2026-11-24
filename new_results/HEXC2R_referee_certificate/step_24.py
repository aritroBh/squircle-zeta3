
step01 = '''"""step_01.py -- referee: validate MY word algebra (reflib) against cmzv6 exactly.
Checks: word lists, shuffle, reg, to/from_series, ds rows, distribution rows,
full relation row SETS at weights 1..4 (primitive integer tuples) -- exact set equality,
plus family counts."""
import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "../HEXC1_hexagon_certificate")
import cmzv6, reflib as R
from fractions import Fraction as Fr
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

# 1. word lists identical
for n in range(0, 5):
    assert R.words(n) == cmzv6.words(n), n
    assert R.words(n, False) == cmzv6.words(n, False), n
log("word lists identical for n=0..4 (counts %s)" % [len(R.words(n)) for n in range(5)]))

# 2. shuffle / reg / series on samples
import random
random.seed(12345)
for trial in range(300):
    n = random.randint(1, 4); m = random.randint(1, 4)
    u = tuple(random.choice(R.LETTERS) for _ in range(n))
    v = tuple(random.choice(R.LETTERS) for _ in range(m))
    assert R.sh(u, v) == {w: Fr(c) for w, c in cmzv6.sh(u, v).items()}, (u, v)
    if u[-1] != -1:
        assert R.to_series(u) == cmzv6.to_series(u)
        assert R.from_series(R.to_series(u)) == u
for w in R.words(3, False)[:]:
    assert R.reg(w) == {x: Fr(c) for x, c in cmzv6.reg(w).items()}, w
log("shuffle/to_series/from_series/reg agree on 300 random pairs + all reg(w), |w|<=3")

# 3. full relation sets: exact set equality at every weight
lower_mine, lower_theirs = {}, {}
for n in (1, 2, 3, 4):
    rows_t, counts_t = cmzv6.make_relations(n, lower_theirs)
    rows_m, counts_m = R.make_relations(n, lower_mine)
    lower_theirs[n] = rows_t; lower_mine[n] = rows_m
    st, sm = set(rows_t), set(rows_m)
    log(f"weight {n}: their rows {len(rows_t)} {counts_t}")
    log(f"          my    rows {len(rows_m)} {counts_m}")
    assert st == sm, f"ROW SETS DIFFER at weight {n}: +{len(st-sm)} -{len(sm-st)}"
    log(f"          SETS EQUAL ({len(st)} rows)")
log("STEP 01 PASS: my relation generators reproduce cmzv6's row sets exactly at weights 1-4")
with open("progress.log", "a") as f:
    f.write("step_01: my word algebra + relation generators match cmzv6 exactly (weights 1-4, set equality)\\n")
'''
with open("steps/step_01.py", "w") as f:
    f.write(step01)
import subprocess
r = subprocess.run(["nice", "-n", "10", "<project>/.venv/bin/python", "steps/step_01.py"], capture_output=True, text=True, timeout=1500)
print(r.stdout)
print(r.stderr[-2000:] if r.returncode else "OK")
