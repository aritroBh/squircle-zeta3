import time, sys
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/HEXC1_hexagon_certificate')
sys.path.insert(0, '<project>/research/08_theory/HEXC2R_referee_certificate')
import cmzv6, reflib as R
from fractions import Fraction as Fr
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

# word lists + algebra spot checks
for n in range(0, 5):
    assert R.words(n) == cmzv6.words(n)
    assert R.words(n, False) == cmzv6.words(n, False)
import random
random.seed(777)
for _ in range(200):
    n = random.randint(1, 4); m = random.randint(1, 4)
    u = tuple(random.choice(R.LETTERS) for _ in range(n))
    v = tuple(random.choice(R.LETTERS) for _ in range(m))
    assert R.sh(u, v) == {w: Fr(c) for w, c in cmzv6.sh(u, v).items()}
log("word lists n<=4 identical; shuffle agrees on 200 random pairs")

# full relation-row SET equality at weights 1..4
lower_m, lower_t = {}, {}
for n in (1, 2, 3, 4):
    rows_t, counts_t = cmzv6.make_relations(n, lower_t)
    rows_m, counts_m = R.make_relations(n, lower_m)
    lower_t[n] = rows_t; lower_m[n] = rows_m
    st, sm = set(rows_t), set(rows_m)
    log(f"weight {n}: cmzv6 {len(rows_t)} rows {counts_t}; reflib {len(rows_m)} rows; SETS EQUAL: {st == sm}")
    assert st == sm
log("STEP PASS: reflib relation generators == cmzv6 row sets exactly at weights 1-4")
with open('progress.log','a') as f:
    f.write("referee step(relsets): reflib make_relations == cmzv6 make_relations as sets, weights 1-4 (4/57/682/7869 rows)\n")
