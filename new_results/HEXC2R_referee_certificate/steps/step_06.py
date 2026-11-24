"""step_06.py -- referee TASK 2 (conversion): my own conversion code (reflib) vs hexc2_common.
(a) EXACT comparison (rational vectors) of all 38 HPL-word conversions and 9 Li-atom
    conversions against hexc2_common (same frames);
(b) NUMERIC audit at 70 dps: my conversion + MY CMZV evaluator vs MY HPL evaluator and
    vs mp.polylog directly;
(c) elementary atoms: my atoms vs their values at 65 dps;
(d) demonstrate the composition order matters: with the WRONG order (outer letters on the
    a-side) the numeric check FAILS -- the audit really discriminates the convention.
"""
import sys, json, re, time
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "../HEXC2_level6_certificate")
sys.path.insert(0, "../HEXC1_hexagon_certificate")
import mpmath as mp
from fractions import Fraction as Fr
import reflib as R
import hexc2_common as HC
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

# the 38 words
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
def load_fit(name):
    D = json.load(open(G1 + "/" + name)); out = []
    for key, val in D:
        m = re.match(r"c\^(-?\d+) H\[(.*)\]", key)
        if key == "1": e, w = 0, ()
        elif key == "c": e, w = 1, ()
        else:
            e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(",")) if m.group(2).strip() else ()
        out.append((e, w))
    return out
WORDS = sorted({w for name in ["hpl_fit_k3_(0,0,0).json", "hpl_fit_k3_(1,0,0).json"]
                for e, w in load_fit(name) if w}, key=lambda x: (len(x), x))

# (a) exact comparison, all 38 words
mismatch = 0
for w in WORDS:
    mine = R.my_hpl_to_cmzv(w)
    theirs = HC.hpl_to_cmzv(w)
    if mine != theirs:
        mismatch += 1
        if mismatch <= 3: print("MISMATCH", w)
log(f"(a) all {len(WORDS)} HPL-word conversions EXACTLY equal to hexc2_common: "
    f"{len(WORDS)-mismatch}/{len(WORDS)}")
assert mismatch == 0

# (a2) exact comparison, 9 Li atoms (same rep=0 frames)
LI = [(2, "1/3"), (3, "1/3"), (3, "-1/2"), (4, "1/3"), (4, "-1/3"), (4, "2/3"),
      (4, "1/2"), (4, "-1/2"), (4, "-2")]
mismatch = 0
for n, z in LI:
    frame = HC.cross_table()[Fr(z)][0]   # same frame as hexc2 for exact comparison
    frame = tuple(None if p is HC.INF else p for p in frame)
    mine = R.my_li_atom(n, Fr(z), frame=frame)
    theirs = HC.li_atom(n, Fr(z))
    if mine != theirs:
        mismatch += 1
        print("MISMATCH Li", n, z)
log(f"(a2) all {len(LI)} Li-atom conversions EXACTLY equal: {len(LI)-mismatch}/{len(LI)}")
assert mismatch == 0

# (b) numeric audit at 70 dps with MY evaluator
mp.mp.dps = 70
worst_hpl = 0; wword = None
for w in WORDS:
    v_conv = R.cmzv_combo(R.my_hpl_to_cmzv(w))
    v_hpl = R.hpl_value(w, mp.mpf(1) / 2)
    err = abs(v_conv - v_hpl)
    if err > worst_hpl: worst_hpl, wword = err, w
log(f"(b) {len(WORDS)} HPL conversions, my evaluator, worst |conv - H_w(1/2)| = "
    f"{mp.nstr(worst_hpl, 3)} (word {wword})")
worst_li = 0; lword = None
for n, z in LI:
    v_conv = R.cmzv_combo(R.my_li_atom(n, Fr(z)))
    v_li = mp.polylog(n, mp.mpf(z))
    err = abs(v_conv - v_li)
    if err > worst_li: worst_li, lword = err, (n, z)
log(f"    {len(LI)} Li atoms, worst |conv - polylog| = {mp.nstr(worst_li, 3)} (Li_{lword[0]}({lword[1]}))")
assert worst_hpl < mp.mpf("1e-62") and worst_li < mp.mpf("1e-62")

# (c) elementary atoms at 65 dps
mp.mp.dps = 65
checks = [
    ("log 2", R.my_atom_log2(), mp.log(2)),
    ("log 3", R.my_atom_log3(), mp.log(3)),
    ("log(3/2)", R.my_atom_log32(), mp.log(mp.mpf(3)/2)),
    ("i pi", R.my_atom_ipi(), mp.pi * mp.j),
    ("zeta(2)", R.my_atom_zeta2(), mp.zeta(2)),
    ("zeta(3)", R.my_atom_zeta3(), mp.zeta(3)),
    ("(i pi)^2 = -pi^2", R.my_atom_pi2(), -mp.pi**2),
    ("(i pi)^4 = pi^4", R.my_atom_pi4(), mp.pi**4),
]
worst_at = 0
for name, vec, v in checks:
    err = abs(R.cmzv_combo(vec) - v)
    worst_at = max(worst_at, err)
    log(f"(c) atom {name:16s}: err {mp.nstr(err, 3)}")
assert worst_at < mp.mpf("1e-58")

# (d) wrong composition order fails the audit
def seg_value_wrong(w, c, e):
    w = tuple(w)
    if c == e: return {(): Fr(1)} if not w else {}
    out = {}; n = len(w)
    for k in range(n + 1):
        left = R.regularize({R._shift6(w[k:], e): Fr(1)})           # WRONG: suffix on the e-side
        right = R.regularize({R._shift6(tuple(reversed(w[:k])), c): Fr(1)})
        out = R.vadd((out, 1), (R.mul(left, right), Fr((-1) ** k)))
    return R.regularize(out)
def hpl_wrong(w, mi=0):
    w = tuple(w)
    formvecs = [R.hpl_form_pullback(mi, sig) for sig in w]
    e0, eh, e1, em, einf = R.MAPS6[mi]
    out = {}
    from itertools import product as P
    for choices in P(*[tuple(fv.items()) for fv in formvecs]):
        ww = tuple(x for cw, a in choices for x in cw)
        coef = Fr(1)
        for cw, a in choices: coef *= a
        out = R.vadd((out, 1), (seg_value_wrong(ww, e0, eh), coef))
    return out
mp.mp.dps = 50
bad = 0; worst_wrong = 0
for w in WORDS:
    err = abs(R.cmzv_combo(hpl_wrong(w)) - R.hpl_value(w, mp.mpf(1)/2))
    worst_wrong = max(worst_wrong, err)
    if err > mp.mpf("1e-30"): bad += 1
log(f"(d) WRONG composition order: {bad}/{len(WORDS)} words off by >1e-30 (worst {mp.nstr(worst_wrong,3)}) "
    "-- audit discriminates the order")
assert bad > 10
log("STEP 06 PASS")
with open("progress.log", "a") as f:
    f.write("step_06: my conversion == hexc2 EXACTLY (38 HPL + 9 Li); numeric 70dps worst %.1e/%.1e; atoms %.1e; wrong order fails (%d/38)\n"
            % (worst_hpl, worst_li, worst_at, bad))
