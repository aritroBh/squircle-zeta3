import time, sys, json, re
sys.dont_write_bytecode = True
base = '<project>/research/08_theory'
sys.path.insert(0, base + '/COR13C_third_corner_certificate')
sys.path.insert(0, base + '/HEXC1_hexagon_certificate')
import mpmath as mp
from fractions import Fraction as Fr
import reflib13 as R13
import reflib as R
import cor13c_common as HC
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

# the 38 HPL words from the G1 fits
WORDS = sorted({w for name in ["hpl_fit_k3_(0,0,0).json", "hpl_fit_k3_(1,0,0).json"]
                for e, w, v in R13.fit_terms(name) if w}, key=lambda x: (len(x), x))
log(f"{len(WORDS)} HPL words from the fits (weights {min(len(w) for w in WORDS)}..{max(len(w) for w in WORDS)})")

# (a) exact comparison, my conversion vs theirs (task map mi=0)
mismatch = 0
for w in WORDS:
    mine = R13.my13_hpl_to_cmzv(w)
    theirs = {w2: Fr(c) for w2, c in HC.hpl_to_cmzv(w).items()}
    if mine != theirs:
        mismatch += 1
        if mismatch <= 3: print("MISMATCH", w)
log(f"(a) HPL-word conversions EXACTLY equal to cor13c_common: {len(WORDS)-mismatch}/{len(WORDS)}")
assert mismatch == 0

# (a2) exact comparison, 9 Li atoms (same frames: HC rep=0 vs my frame matching HC's choice)
LI = [(2, "1/3"), (2, "-1/3"), (3, "1/3"), (3, "-1/2"), (4, "1/3"), (4, "-1/3"),
      (4, "2/3"), (4, "1/2"), (4, "-2")]
mismatch = 0
for n, z in LI:
    frame = HC.cross_table()[Fr(z)][0]
    frame = tuple(None if p is HC.INF else p for p in frame)
    mine = R13.my13_li_atom(n, Fr(z), frame=frame)
    theirs = {w2: Fr(c) for w2, c in HC.li_atom(n, Fr(z)).items()}
    if mine != theirs:
        mismatch += 1
        print("MISMATCH Li", n, z)
log(f"(a2) Li-atom conversions EXACTLY equal: {len(LI)-mismatch}/{len(LI)}")
assert mismatch == 0

# (a3) targets exact equality
T_my = {2: R13.my_target2(), 3: R13.my_target3(), 4: R13.my_target4()}
T_hc = {2: {w: Fr(c) for w, c in HC.target2().items()},
        3: {w: Fr(c) for w, c in HC.target3().items()},
        4: {w: Fr(c) for w, c in HC.target4().items()}}
for n in (2, 3, 4):
    same = T_my[n] == T_hc[n]
    log(f"(a3) my target{n} == cor13c target{n}: {same} (supports {len(T_my[n])}/{len(T_hc[n])})")
    assert same

# (b) numeric audit at 70 dps with MY evaluator
mp.mp.dps = 70
worst_hpl = 0; wword = None
for w in WORDS:
    v_conv = R.cmzv_combo(R13.my13_hpl_to_cmzv(w))
    v_hpl = R.hpl_value(w, mp.mpf(1)/3)
    err = abs(v_conv - v_hpl)
    if err > worst_hpl: worst_hpl, wword = err, w
log(f"(b) {len(WORDS)} HPL conversions, worst |conv - H_w(1/3)| = {mp.nstr(worst_hpl,3)} (word {wword})")
worst_li = 0; latom = None
for n, z in LI:
    v_conv = R.cmzv_combo(R13.my13_li_atom(n, Fr(z)))
    err = abs(v_conv - mp.polylog(n, mp.mpf(z)))
    if err > worst_li: worst_li, latom = err, (n, z)
log(f"    {len(LI)} Li atoms, worst |conv - polylog| = {mp.nstr(worst_li,3)} (Li_{latom[0]}({latom[1]}))")
assert worst_hpl < mp.mpf('1e-62') and worst_li < mp.mpf('1e-62')
with open('progress.log','a') as f:
    f.write("referee step(conversion): my conversion == cor13c EXACTLY (38 HPL + 9 Li + targets); numeric 70dps worst %.1e/%.1e\n" % (worst_hpl, worst_li))
log("STEP PASS (exact equality + 70 dps audit part 1)")
