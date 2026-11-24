import time, sys, random
sys.dont_write_bytecode = True
base = '<project>/research/08_theory'
sys.path.insert(0, base + '/COR13C_third_corner_certificate')
sys.path.insert(0, base + '/HEXC2R_referee_certificate')
import mpmath as mp
mp.mp.dps = 65
import fast_cmzv as FC
import numeric_cmzv as NC
import reflib as R
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

FC.setup(62); FC.finish_setup()
NC.set_dps(62)

random.seed(20250)
words_tested = set()
worst_fc = worst_nc = 0
ntest = 0
while ntest < 150:
    n = random.randint(1, 4)
    w = tuple(random.choice([-1,0,1,2,3,4,5]) for _ in range(n))
    if w and (w[0] == 0 or w[-1] == -1): continue
    if w in words_tested: continue
    words_tested.add(w)
    a = R.cmzv_value(w)
    b = FC.cmzv_value(w); b = mp.mpc(b.real.mid(), b.imag.mid())
    c = NC.cmzv_value(w)
    worst_fc = max(worst_fc, abs(a - b))
    worst_nc = max(worst_nc, abs(a - c))
    ntest += 1
log(f"150 random admissible words (weights 1-4): |reflib - fast_cmzv| worst {mp.nstr(worst_fc,3)}, |reflib - numeric_cmzv| worst {mp.nstr(worst_nc,3)}")
assert worst_fc < mp.mpf('1e-58') and worst_nc < mp.mpf('1e-58')

# known constants as words: zeta(2)=(-1,0), zeta(3)=(-1,-1,0), Li_2(zeta)=(1,1)? value of word (j..): check a few against polylogs
Z = mp.e**(mp.pi*mp.j/3)
checks = [
    ((-1, 0), mp.zeta(2)),
    ((-1, -1, 0), mp.zeta(3)),
    ((-1, -1, -1, 0), mp.zeta(4)),
    ((0, 0, 1), mp.polylog(3, 1) if False else None),
    ((1,), -mp.log(1 - Z)),
    ((2,), -mp.log(1 - Z**2)),
    ((3,), -mp.log(2)),
    ((0, 1), mp.polylog(2, 1)),  # word (0,1): int dt/t then dt/(1-t) = Li_2(1)? careful: letter 0 = zeta^0 dt/(1-t)=dt/(1-t)? no: letter j=0 is dt/(1-t); so (0,1) is Li_2(zeta)
]
for w, v in checks:
    if v is None: continue
    err = abs(R.cmzv_value(w) - v)
    log(f"  word {w}: |reflib - known| = {mp.nstr(err,3)}")
    assert err < mp.mpf('1e-58')
log("evaluator cross-validation PASS")
with open('progress.log','a') as f:
    f.write("referee step(eval xcheck): reflib evaluator vs fast_cmzv/numeric_cmzv on 150 words, worst %.1e/%.1e; constants OK\n" % (worst_fc, worst_nc))
