import time, sys
sys.dont_write_bytecode = True
base = '<project>/research/08_theory'
sys.path.insert(0, base + '/HEXC2R_referee_certificate')
import mpmath as mp
from fractions import Fraction as Fr
import reflib13 as R13
import reflib as R
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

WORDS = sorted({w for name in ["hpl_fit_k3_(0,0,0).json", "hpl_fit_k3_(1,0,0).json"]
                for e, w, v in R13.fit_terms(name) if w}, key=lambda x: (len(x), x))
LI = [(2, "1/3"), (2, "-1/3"), (3, "1/3"), (3, "-1/2"), (4, "1/3"), (4, "-1/3"),
      (4, "2/3"), (4, "1/2"), (4, "-2")]

# (b) numeric audit at 70 dps, MY conversion + MY evaluator vs MY hpl_value / polylog
mp.mp.dps = 70
worst_hpl = 0; wword = None
for w in WORDS:
    err = abs(R.cmzv_combo(R13.my13_hpl_to_cmzv(w)) - R.hpl_value(w, mp.mpf(1)/3))
    if err > worst_hpl: worst_hpl, wword = err, w
log(f"(b) {len(WORDS)} HPL conversions at 1/3, worst |conv - H_w(1/3)| = {mp.nstr(worst_hpl,3)} (word {wword})")
worst_li = 0; latom = None
for n, z in LI:
    err = abs(R.cmzv_combo(R13.my13_li_atom(n, Fr(z))) - mp.polylog(n, mp.mpf(z)))
    if err > worst_li: worst_li, latom = err, (n, z)
log(f"    {len(LI)} Li atoms, worst |conv - polylog| = {mp.nstr(worst_li,3)} (Li_{latom[0]}({latom[1]}))")
assert worst_hpl < mp.mpf('1e-62') and worst_li < mp.mpf('1e-62')

# (c) elementary atoms at 70 dps
checks = [
    ("log 2", R13.a_log2(), mp.log(2)),
    ("log 3", R13.a_log3(), mp.log(3)),
    ("i pi", R13.a_ipi(), mp.pi * mp.j),
    ("zeta(2)", R13.a_zeta2(), mp.zeta(2)),
    ("zeta(3)", R13.a_zeta3(), mp.zeta(3)),
    ("(i pi)^2 = -pi^2", R13.a_pi2(), -mp.pi**2),
    ("(i pi)^4 = pi^4", R13.a_pi4(), mp.pi**4),
]
worst_at = 0
for name, vec, v in checks:
    err = abs(R.cmzv_combo(vec) - v)
    worst_at = max(worst_at, err)
    log(f"    atom {name:16s}: |err| = {mp.nstr(err,3)}")
assert worst_at < mp.mpf('1e-62')

# (d) WRONG composition order fails: suffix on the e-side
def seg_wrong(w, c, e):
    w = tuple(w)
    if c == e: return {(): Fr(1)} if not w else {}
    out = {}; n = len(w)
    for k in range(n + 1):
        left = R.regularize({R._shift6(w[k:], e): Fr(1)})
        right = R.regularize({R._shift6(tuple(reversed(w[:k])), c): Fr(1)})
        out = R.vadd((out, 1), (R.mul(left, right), Fr((-1) ** k)))
    return R.regularize(out)
def hpl_wrong(w, mi=0):
    from itertools import product as P
    w = tuple(w)
    formvecs = [R13.hpl_form_pullback13(mi, sig) for sig in w]
    e0, eh, e1, em, einf = R13.maps13()[mi]
    out = {}
    for choices in P(*[tuple(fv.items()) for fv in formvecs]):
        ww = tuple(x for cw, a in choices for x in cw)
        coef = Fr(1)
        for cw, a in choices: coef *= a
        out = R.vadd((out, 1), (seg_wrong(ww, e0, eh), coef))
    return out
bad = 0; worst_wrong = 0
for w in WORDS:
    err = abs(R.cmzv_combo(hpl_wrong(w)) - R.hpl_value(w, mp.mpf(1)/3))
    worst_wrong = max(worst_wrong, err)
    if err > mp.mpf('1e-30'): bad += 1
log(f"(d) WRONG composition order: {bad}/{len(WORDS)} words off by >1e-30 (worst {mp.nstr(worst_wrong,3)}) -- audit discriminates the order")
assert bad > 10
with open('progress.log','a') as f:
    f.write("referee step(conversion audit): 70dps HPL worst %.1e, Li worst %.1e, atoms worst %.1e; wrong order fails %d/38\n" % (worst_hpl, worst_li, worst_at, bad))
log("STEP PASS")
