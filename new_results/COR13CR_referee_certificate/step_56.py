import time, sys, json
sys.dont_write_bytecode = True
base = '<project>/research/08_theory'
sys.path.insert(0, base + '/HEXC1_hexagon_certificate')
sys.path.insert(0, base + '/HEXC2R_referee_certificate')
sys.path.insert(0, base + '/COR13C_third_corner_certificate')
import mpmath as mp
from fractions import Fraction as Fr
import cmzv6
import reflib as R
import reflib13 as R13
from cor13c_common import target3
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

# (1) TAMPER TEST at weight 3
cert23 = json.load(open(base + '/COR13C_third_corner_certificate/cert13_w23.json'))
lower = {}
for n in (1, 2, 3):
    rows, _ = cmzv6.make_relations(n, lower); lower[n] = rows
C = cert23["3"]
prow, y = C["prow"], [Fr(int(a), int(b)) for a, b in C["y"]]
T = {w: Fr(c) for w, c in target3().items()}
cols = list(cmzv6.words(3))
def resid_count(sel, yv, T, cols):
    ids = {w: j for j, w in enumerate(cols)}
    acc = [Fr(0)]*len(cols)
    for row, yi in zip(sel, yv):
        if yi == 0: continue
        for w, c in row:
            acc[ids[w]] += yi*Fr(c)
    return sum(1 for w, j in ids.items() if acc[j] != T.get(w, Fr(0)))
sel = [lower[3][i] for i in prow]
assert resid_count(sel, y, T, cols) == 0
y_bad = list(y); y_bad[0] += 1
bad = resid_count(sel, y_bad, T, cols)
log(f"(1) tamper y[0]+=1: residual words = {bad} (must be > 0)")
assert bad > 0
y_bad2 = list(y); y_bad2[100] = 2*y_bad2[100]
bad2 = resid_count(sel, y_bad2, T, cols)
log(f"    tamper y[100]*=2: residual words = {bad2} (must be > 0)")
assert bad2 > 0

# (2) MAP INDEPENDENCE: targets under other maps vanish numerically
mp.mp.dps = 60
for mi in (3, 7, 11):
    v2 = R.cmzv_combo(R13.my_target2(mi))
    log(f"(2) map mi={mi} (exps {R13.maps13()[mi]}): |target2 value| = {mp.nstr(abs(v2),3)}")
    assert abs(v2) < mp.mpf('1e-50')
v3 = R.cmzv_combo(R13.my_target3(5))
log(f"    map mi=5 (exps {R13.maps13()[5]}): |target3 value| = {mp.nstr(abs(v3),3)}")
assert abs(v3) < mp.mpf('1e-50')
log("STEP PASS: tamper test discriminates; targets vanish under other maps too")
with open('progress.log','a') as f:
    f.write("referee step(tamper+maps): tampered certs FAIL (%d,%d words); target2 vanishes under maps 3,7,11; target3 under map 5\n" % (bad, bad2))
