code = '''"""verify_hexc2.py -- HEXC2 certificate verification. Runtime target < 5 min.

PART A (the proof, exact rational arithmetic):
  For each weight n in {2,3,4}: rebuild the proved level-6 relation rows with
  cmzv6.make_relations (regularized double shuffle [Ihara-Kaneko-Zagier], Zhao's
  regularized distribution d in {2,3,6}, weight-1 branch identities, shuffle lifts),
  rebuild the target vector T_n (ID_n: hexagon-corner side minus closed form, converted
  to level-6 CMZV words via the Mobius map f), load the certificate (row indices prow,
  rational multipliers y) from cert_w23.json / cert_w4.json, and check
      sum_i y_i rows[prow_i] - T_n == 0   exactly, on every admissible word.
  Then T_n is an exact rational combination of proved relations, hence evaluates to 0.

PART B (numeric audit of the conversion code, 62 dps; not part of the proof):
  Mobius map images, all 38 HPL-word conversions vs H_series, all 9 Li atoms vs polylog,
  elementary atoms, the three target values, the Li4(1/9) distribution input,
  and a cross-check of the two independent CMZV evaluators.
"""
import os, sys, json, time, re
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
THEORY = "<project>/research/08_theory"
sys.path.insert(0, THEORY + "/HEXC1_hexagon_certificate")
sys.path.insert(0, THEORY + "/COR4_M3alpha_classical")
from fractions import Fraction as Fr
T0 = time.monotonic()
def log(s):
    print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

PART_B = ("--fast" not in sys.argv)

# ---------------- Part A ----------------
log("PART A: exact certificate verification")
import cmzv6
from hexc2_common import target2, target3, target4
from flint import fmpq_mat

cert23 = json.load(open(HERE + "/cert_w23.json"))
cert4 = json.load(open(HERE + "/cert_w4.json"))

def exact_verify(n, rows, cols, prow, y, T):
    sel = [rows[i] for i in prow]
    ids = {w: i for i, w in enumerate(cols)}
    assert set(T) <= set(cols), "target word outside admissible words"
    r, m = len(sel), len(cols)
    S = fmpq_mat(r, m)
    for i, row in enumerate(sel):
        for w, c in row:
            S[i, ids[w]] = c
    Y = fmpq_mat(r, 1)
    for i, c in enumerate(y):
        Y[i, 0] = fmpq_mat(1, 1, [str(c)])[0, 0]
    res = Y.transpose() * S
    bad = 0
    for w, j in ids.items():
        if Fr(str(res[0, j])) != T.get(w, Fr(0)):
            bad += 1
            if bad < 3: print("   MISMATCH at", w)
    return bad

lower = {}
RELS = {}
for n in (1, 2, 3, 4):
    rows, counts = cmzv6.make_relations(n, lower)
    lower[n] = rows
    RELS[n] = rows
    log(f"  built relations weight {n}: {len(rows)} rows, families {counts}")
    if n < 2: continue
    T = {2: target2, 3: target3, 4: target4}[n]()
    cols = cmzv6.words(n)
    if n < 4:
        C = cert23[str(n)]
        prow, y = C["prow"], [Fr(int(a), int(b)) for a, b in C["y"]]
    else:
        prow, y = cert4["prow"], [Fr(s) for s in cert4["y"]]
    bad = exact_verify(n, rows, cols, prow, y, T)
    status = "PASS" if bad == 0 else f"FAIL ({bad} mismatches)"
    log(f"  weight {n}: certificate rows={len(prow)} target_support={len(T)} "
        f"residual_words={bad}  => {status}")
    assert bad == 0
log("PART A PASSED: ID2, ID3, ID4 are exact rational combinations of proved relations")

if not PART_B:
    log("(--fast: skipping numeric audit)")
    raise SystemExit(0)

# ---------------- Part B ----------------
log("PART B: numeric audit at 62 dps (validates conversion code; the proof is Part A)")
import mpmath as mp
import fast_cmzv as FC
import numeric_cmzv as NC
from hexc2_common import (hpl_to_cmzv, li_atom, atom_log2, atom_log3, atom_ipi,
                          atom_zeta2, atom_zeta3, atom_pi2, atom_pi4, MAPS)
from hpl_series import H_series
mp.mp.dps = 65
FC.setup(62); FC.finish_setup()

def combo(d):
    tot = mp.mpc(0)
    for w, c in d.items():
        v = FC.cmzv_value(w)
        tot += mp.mpf(c.numerator) / c.denominator * mp.mpc(v.real.mid(), v.imag.mid())
    return tot

ok = True
def check(name, resid, tol):
    global ok
    good = abs(resid) < tol
    ok = ok and good
    log(f"  [{'PASS' if good else 'FAIL'}] {name}: |resid| = {mp.nstr(abs(resid), 3)}")

# B0. The 12 Mobius maps send {0,1/2,1,-1,inf} to the claimed sixth roots (numeric).
zeta = mp.e ** (mp.pi * mp.j / 3)
# task map f(z) = (zeta^2 z + 1)/(zeta^4 z + 1) has exps (0,1,2,5,4)
f = lambda z: (zeta**2 * z + 1) / (zeta**4 * z + 1)
imgs = [f(0), f(mp.mpf(1)/2), f(1), f(-1)]
err = max(abs(imgs[i] - zeta**e) for i, e in enumerate((0, 1, 2, 5)))
check("task Mobius map f on {0,1/2,1,-1}", err, mp.mpf("1e-40"))
check("task Mobius map f(inf) = zeta^4", abs(zeta**2/zeta**4 - zeta**4), mp.mpf("1e-40"))

# B1. every HPL-word conversion vs H_series
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
def load_fit(name):
    D = json.load(open(G1 + "/" + name)); out = []
    for key, val in D:
        m = re.match(r"c\\^(-?\\d+) H\\[(.*)\\]", key)
        if key == "1": e, w = 0, ()
        elif key == "c": e, w = 1, ()
        else:
            e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(",")) if m.group(2).strip() else ()
        out.append((e, w))
    return out
WORDS = sorted({w for name in ["hpl_fit_k3_(0,0,0).json", "hpl_fit_k3_(1,0,0).json"]
                for e, w in load_fit(name) if w}, key=lambda x: (len(x), x))
worst = 0; wword = None
for w in WORDS:
    err = abs(combo(hpl_to_cmzv(w)) - H_series(w, mp.mpf(1) / 2))
    if err > worst: worst, wword = err, w
check(f"all {len(WORDS)} HPL-word conversions (worst: H{list(wword)})", worst, mp.mpf("1e-58"))

# B2. every Li atom vs polylog
LI = [(2, "1/3"), (3, "1/3"), (3, "-1/2"), (4, "1/3"), (4, "-1/3"), (4, "2/3"),
      (4, "1/2"), (4, "-1/2"), (4, "-2")]
worst = 0; latom = None
for n, z in LI:
    err = abs(combo(li_atom(n, Fr(z))) - mp.polylog(n, mp.mpf(z)))
    if err > worst: worst, latom = err, (n, z)
check(f"all {len(LI)} Li-atom conversions (worst: Li_{latom[0]}({latom[1]}))", worst, mp.mpf("1e-58"))

# B3. elementary atoms
check("zeta(2) atom", abs(combo(atom_zeta2()) - mp.zeta(2)), mp.mpf("1e-58"))
check("zeta(3) atom", abs(combo(atom_zeta3()) - mp.zeta(3)), mp.mpf("1e-58"))
check("log2 atom", abs(combo(atom_log2()) - mp.log(2)), mp.mpf("1e-58"))
check("log3 atom", abs(combo(atom_log3()) - mp.log(3)), mp.mpf("1e-58"))
check("i*pi atom", abs(combo(atom_ipi()) - mp.pi * mp.j), mp.mpf("1e-58"))
check("(i*pi)^2 atom = -pi^2", abs(combo(atom_pi2()) + mp.pi**2), mp.mpf("1e-58"))
check("(i*pi)^4 atom = pi^4", abs(combo(atom_pi4()) - mp.pi**4), mp.mpf("1e-58"))

# B4. the three targets
for name, T in [("ID2", target2), ("ID3", target3), ("ID4", target4)]:
    v = combo(T())
    check(f"target {name} value (support {len(T())})", v, mp.mpf("1e-55"))

# B5. input theorem: Li4(1/9) = 8 (Li4(1/3) + Li4(-1/3))  (absolutely convergent series)
lhs = mp.polylog(4, mp.mpf(1) / 9) - 8 * (mp.polylog(4, mp.mpf(1) / 3) + mp.polylog(4, -mp.mpf(1) / 3))
check("Li4(1/9) distribution (input)", lhs, mp.mpf("1e-58"))

# B6. cross-check of the two independent CMZV evaluators
NC.set_dps(62)
import random
random.seed(7)
worst = 0
for _ in range(10):
    w = tuple(random.choice([-1, 0, 1, 2, 3, 4, 5]) for _ in range(4))
    if w[0] == 0 or w[-1] == -1: continue
    a = FC.cmzv_value(w); b = NC.cmzv_value(w)
    err = abs(mp.mpc(a.real.mid(), a.imag.mid()) - b)
    worst = max(worst, err)
check("two independent CMZV evaluators (10 random words)", worst, mp.mpf("1e-55"))

log("PART B " + ("PASSED" if ok else "FAILED"))
log("VERDICT: " + ("HEXC2 CERTIFICATE VALID — ID2, ID3, ID4 PROVED" if ok else "FAILURE"))
assert ok
'''
open('verify_hexc2.py','w').write(code)
print('written verify_hexc2.py', len(code))
