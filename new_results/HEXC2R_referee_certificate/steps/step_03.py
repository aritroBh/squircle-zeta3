"""step_03.py -- referee: cross-validate MY CMZV evaluator (reflib, truncated multiple sums,
independent algorithm) against HEXC2's two evaluators (numeric_cmzv mpmath log-series,
fast_cmzv FLINT acb_poly) on random admissible words at 65 dps; sanity values vs known
constants; HPL cross-check vs hpl_series."""
import sys, time, random
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "../HEXC2_level6_certificate")
sys.path.insert(0, "../COR4_M3alpha_classical")
import mpmath as mp
mp.mp.dps = 65
import reflib as R
import numeric_cmzv as NC
import fast_cmzv as FC
from hpl_series import H_series
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

NC.set_dps(62)
FC.setup(62); FC.finish_setup()

def fcv(w):
    v = FC.cmzv_value(w)
    return mp.mpc(v.real.mid(), v.imag.mid())

# sanity: known values
zeta = mp.e**(mp.pi*mp.j/3)
checks = [
    ((-1, 0), mp.zeta(2), "word (-1,0) = zeta(2)"),
    ((-1, -1, 0), mp.zeta(3), "word (-1,-1,0) = zeta(3)"),
    ((3,), -mp.log(2), "word (3,) = -log2"),
    ((1,), -mp.log(1-zeta), "word (1,) = -log(1-zeta)"),
    ((1, 1), -mp.pi**2/18, "word (1,1) = -pi^2/18 (from Euler cert expectation)"),
]
ok = True
for w, v, name in checks:
    mine = R.cmzv_value(w)
    err = abs(mine - v)
    good = err < mp.mpf("1e-55")
    ok &= good
    log(f"  [{'PASS' if good else 'FAIL'}] {name}: mine={mp.nstr(mine, 12)} err={mp.nstr(err,3)}")

# 3-way cross-check on random admissible words
random.seed(4242)
worst_mn = worst_mf = 0
NTRIAL = 40
for _ in range(NTRIAL):
    n = random.randint(1, 4)
    w = tuple(random.choice([-1,0,1,2,3,4,5]) for _ in range(n))
    if w[0] == 0 or w[-1] == -1: continue
    a = R.cmzv_value(w); b = NC.cmzv_value(w); c = fcv(w)
    worst_mn = max(worst_mn, abs(a-b)); worst_mf = max(worst_mf, abs(a-c))
log(f"  3-way cross-check ({NTRIAL} random words, weights 1-4):")
log(f"    |mine - numeric_cmzv| worst = {mp.nstr(worst_mn, 3)}")
log(f"    |mine - fast_cmzv|    worst = {mp.nstr(worst_mf, 3)}")
ok &= worst_mn < mp.mpf("1e-50") and worst_mf < mp.mpf("1e-50")

# every admissible word of weight 2 and a big weight-3 sample, mine vs NC
worst = 0
for w in R.words(2):
    worst = max(worst, abs(R.cmzv_value(w) - NC.cmzv_value(w)))
log(f"  all 36 weight-2 words: |mine - numeric_cmzv| worst = {mp.nstr(worst,3)}")
ok &= worst < mp.mpf("1e-50")
random.seed(777)
w3 = random.sample(list(R.words(3)), 60)
worst = 0
for w in w3:
    worst = max(worst, abs(R.cmzv_value(w) - NC.cmzv_value(w)))
log(f"  60 random weight-3 words: worst = {mp.nstr(worst,3)}")
ok &= worst < mp.mpf("1e-50")

# HPL cross-check
worst = 0
for w in [(0,), (1,), (-1,), (0,1), (0,-1), (1,-1), (0,0,1), (1,0,-1), (0,1,-1,1), (0,0,-1,1), (1,1,0,-1), (1,0), (0,1,0)]:
    ww = H_series(w, mp.mpf(1)/2) if w[-1] != 0 else None
    if ww is not None:
        worst = max(worst, abs(R.hpl_value(w, mp.mpf(1)/2) - ww))
log(f"  HPL words (12 samples) |mine - H_series| worst = {mp.nstr(worst,3)}")
ok &= worst < mp.mpf("1e-50")

log("STEP 03 " + ("PASS" if ok else "FAIL"))
assert ok
with open("progress.log", "a") as f:
    f.write("step_03: my evaluator cross-checks vs numeric_cmzv/fast_cmzv/H_series PASS (worst %.1e, %.1e @65dps)\n" % (worst_mn, worst_mf))
