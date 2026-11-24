import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, '.')
import mpmath as mp
import fast_cmzv as FC
import numeric_cmzv as NC
from cor13c_common import target2, target3, target4
t0 = time.monotonic()
FC.setup(62); FC.finish_setup()
mp.mp.dps = 65

def combo(d):
    tot = mp.mpc(0)
    for w, c in d.items():
        v = FC.cmzv_value(w)
        tot += mp.mpf(c.numerator) / c.denominator * mp.mpc(v.real.mid(), v.imag.mid())
    return tot

res = {}
for name, T in [("ID2", target2), ("ID3", target3), ("ID4", target4)]:
    TT = T()
    v = combo(TT)
    res[name] = (len(TT), abs(v))
    print(f"[target {name}] support {len(TT)}, |value| = {mp.nstr(abs(v),3)} ({time.monotonic()-t0:.1f}s)")

# cross-check the two evaluators on 12 random weight-4 words
NC.set_dps(62)
import random
random.seed(13)
worst = 0
for _ in range(12):
    w = tuple(random.choice([-1,0,1,2,3,4,5]) for _ in range(4))
    if w[0] == 0 or w[-1] == -1: continue
    a = FC.cmzv_value(w); b = NC.cmzv_value(w)
    err = abs(mp.mpc(a.real.mid(), a.imag.mid()) - b)
    worst = max(worst, err)
print(f"[evaluators] FC vs NC on 12 random words, worst |diff| = {mp.nstr(worst,3)}")
ok = all(v < mp.mpf('1e-55') for _, v in res.values()) and worst < mp.mpf('1e-55')
print("STEP_25", "PASS" if ok else "FAIL")
with open('progress.log','a') as f:
    f.write(f"step_25: targets at 62dps {[(k, mp.nstr(v,2)) for k,(s,v) in res.items()]}; evaluators {mp.nstr(worst,2)} -> {'PASS' if ok else 'FAIL'}\n")
