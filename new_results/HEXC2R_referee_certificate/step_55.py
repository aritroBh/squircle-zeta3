
step07 = '''"""step_07.py -- referee: rebuild T2/T3/T4 with MY OWN conversion (reflib) and MY OWN
transcription of the COR5 closed forms, exact-compare with hexc2_common's targets, and
re-run the exact certificate check against MY targets. Also resolve the sign of H_{0,-1}
in the printed ID2 (result.md of HEXC2 says -H_{0,-1}; the fit says +H_{0,-1})."""
import sys, json, re, time
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "../HEXC2_level6_certificate")
sys.path.insert(0, "../HEXC1_hexagon_certificate")
import mpmath as mp
from fractions import Fraction as Fr
import reflib as R
import hexc2_common as HC
import cmzv6
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
def my_fit_vec(name, weight):
    """sum cf * (1/2)^e * H_w(1/2) over terms of HPL-weight `weight`, MY conversion."""
    D = json.load(open(G1 + "/" + name))
    out = {}
    for key, val in D:
        m = re.match(r"c\\^(-?\\d+) H\\[(.*)\\]", key)
        if key == "1": e, w = 0, ()
        elif key == "c": e, w = 1, ()
        else:
            e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(",")) if m.group(2).strip() else ()
        if len(w) != weight: continue
        v = {(): Fr(1)} if not w else R.my_hpl_to_cmzv(w)
        out = R.vadd((out, 1), (v, Fr(val) * Fr(1, 2) ** e))
    return out

def prod(*vecs):
    out = {(): Fr(1)}
    for v in vecs: out = R.mul(out, v)
    return out

def Li(n, z): return R.my_li_atom(n, Fr(z))

def my_F4_hpl_vec():
    D = json.load(open(G1 + "/hpl_fit_k3_(0,0,0).json"))
    out = {}
    for key, val in D:
        m = re.match(r"c\\^(-?\\d+) H\\[(.*)\\]", key)
        if not m: continue
        e = int(m.group(1))
        w = tuple(int(z) for z in m.group(2).split(",")) if m.group(2).strip() else ()
        if len(w) != 4: continue
        out = R.vadd((out, 1), (R.my_hpl_to_cmzv(w), Fr(val) * Fr(1, 2) ** e))
    return out

def my_F4_closed():
    L2, L3 = R.my_atom_log2(), R.my_atom_log3()
    Z3, PI2, PI4 = R.my_atom_zeta3(), {w: -c for w, c in R.my_atom_pi2().items()}, R.my_atom_pi4()
    Li21, Li31, Li41 = Li(2, '1/3'), Li(3, '1/3'), Li(4, '1/3')
    Li3m, Li4m = Li(3, '-1/2'), Li(4, '-1/2')
    out = {}
    def T(vec, coef):
        nonlocal out; out = R.vadd((out, 1), (vec, Fr(coef)))
    T(Li(4, -2), -21); T(Li41, Fr(-65, 8)); T(Li(4, '2/3'), Fr(13, 2)); T(Li(4, '1/2'), Fr(137, 6))
    T(Li41, Fr(39, 8)); T(Li(4, '-1/3'), Fr(39, 8))      # 39/64 Li4(1/9) by distribution
    T(prod(Li31, L2), Fr(15, 2)); T(prod(Li3m, L2), Fr(37, 2)); T(prod(Li3m, L3), Fr(-5, 4))
    T(prod(Z3, L2), Fr(119, 16)); T(prod(Z3, L3), Fr(-37, 32))
    T(prod(Li21, Li21), 1); T(prod(Li21, L2, L2), Fr(-23, 4)); T(prod(Li21, L2, L3), Fr(7, 4))
    T(prod(Li21, L3, L3), 1); T(prod(PI2, Li21), Fr(-1, 6))
    T(prod(PI2, L2, L2), Fr(-455, 144)); T(prod(PI2, L2, L3), Fr(9, 8)); T(prod(PI2, L3, L3), Fr(-7, 32))
    T(PI4, Fr(-949, 2160))
    T(prod(L2, L2, L2, L2), Fr(-211, 144)); T(prod(L2, L2, L2, L3), Fr(101, 24))
    T(prod(L2, L2, L3, L3), Fr(-15, 8)); T(prod(L2, L3, L3, L3), Fr(-35, 24))
    T(prod(L3, L3, L3, L3), Fr(113, 192))
    return out

def my_W3_closed():
    L2, L3 = R.my_atom_log2(), R.my_atom_log3()
    Z3, PI2 = R.my_atom_zeta3(), {w: -c for w, c in R.my_atom_pi2().items()}
    Li21, Li31, Li3m = Li(2, '1/3'), Li(3, '1/3'), Li(3, '-1/2')
    out = {}
    def T(vec, coef):
        nonlocal out; out = R.vadd((out, 1), (vec, Fr(coef)))
    T(Z3, Fr(-15, 4)); T(Li31, -3); T(Li3m, -8)
    T(prod(PI2, L2), Fr(-1, 12)); T(prod(Li21, L2), 4); T(prod(Li21, L3), Fr(-3, 2))
    T(prod(L2, L2, L2), Fr(1, 3)); T(prod(L2, L2, L3), Fr(-5, 2)); T(prod(L2, L3, L3), 2)
    T(prod(L3, L3, L3), Fr(-1, 4))
    return out

def my_W2_closed():
    L2, L32 = R.my_atom_log2(), R.my_atom_log32()
    return R.vadd((prod(L2, L2), 1), (Li(2, '1/3'), 1), (prod(L32, L32), 1)) if False else \
        R.vadd((prod(L2, L2), Fr(1, 2)), (Li(2, '1/3'), 1), (prod(L32, L32), Fr(1, 2)))

def my_target2():
    R32 = my_fit_vec("hpl_fit_k3_(0,0,0).json", 2)
    return R.vadd((R32, 1), (R.my_atom_zeta2(), Fr(1, 2)), (my_W2_closed(), -1))

def my_target3():
    R33 = my_fit_vec("hpl_fit_k3_(0,0,0).json", 3)
    Psi1 = my_fit_vec("hpl_fit_k3_(1,0,0).json", 1)
    lhs = R.vadd((R33, 1), (R.mul(R.my_atom_zeta2(), Psi1), 1), (R.my_atom_zeta3(), Fr(-1, 2)))
    return R.vadd((lhs, 1), (my_W3_closed(), -1))

def my_target4():
    return R.vadd((my_F4_hpl_vec(), 1), (my_F4_closed(), -1))

# (1) exact-compare against hexc2_common targets
T_my = {2: my_target2(), 3: my_target3(), 4: my_target4()}
T_hc = {2: HC.target2(), 3: HC.target3(), 4: HC.target4()}
for n in (2, 3, 4):
    same = T_my[n] == T_hc[n]
    log(f"(1) my target{n} == hexc2 target{n}: {same} (support {len(T_my[n])})")
    assert same

# (2) exact certificate check against MY targets
from flint import fmpq_mat
cert23 = json.load(open("../HEXC2_level6_certificate/cert_w23.json"))
cert4 = json.load(open("../HEXC2_level6_certificate/cert_w4.json"))
lower = {}
for n in (1, 2, 3, 4):
    rows, counts = cmzv6.make_relations(n, lower)
    lower[n] = rows
    if n < 2: continue
    T = T_my[n]
    cols = list(cmzv6.words(n))
    ids = {w: j for j, w in enumerate(cols)}
    if n < 4:
        C = cert23[str(n)]
        prow, y = C["prow"], [Fr(int(a), int(b)) for a, b in C["y"]]
    else:
        prow, y = cert4["prow"], [Fr(s) for s in cert4["y"]]
    sel = [rows[i] for i in prow]
    Sm = fmpq_mat(len(sel), len(cols))
    for i, row in enumerate(sel):
        for w, c in row:
            Sm[i, ids[w]] = c
    Y = fmpq_mat(len(sel), 1)
    for i, yi in enumerate(y):
        Y[i, 0] = fmpq_mat(1, 1, [str(yi)])[0, 0]
    res = Y.transpose() * Sm
    bad = sum(1 for w, j in ids.items() if Fr(str(res[0, j])) != T.get(w, Fr(0)))
    log(f"(2) weight {n}: cert vs MY target, residual_words={bad}  {'PASS' if bad == 0 else 'FAIL'}")
    assert bad == 0

# (3) numeric values of my targets at 65 dps
mp.mp.dps = 65
for n in (2, 3, 4):
    v = R.cmzv_combo(T_my[n])
    log(f"(3) |my target{n} value| = {mp.nstr(abs(v), 3)}")
    assert abs(v) < mp.mpf("1e-55")

# (4) the printed ID2 sign: numeric resolution
mp.mp.dps = 50
H01 = R.hpl_value((0, 1), mp.mpf(1)/2)
H0m = R.hpl_value((0, -1), mp.mpf(1)/2)
rhs = mp.log(2)**2/2 + mp.polylog(2, mp.mpf(1)/3) + mp.log(mp.mpf(3)/2)**2/2
vA = -H01 + H0m + mp.zeta(2)/2     # fit coefficients: -H[0,1] + H[0,-1]
vB = -H01 - H0m + mp.zeta(2)/2     # as printed in HEXC2 result.md ID2
log(f"(4) -H01 + H0-1 + zeta(2)/2 = {mp.nstr(vA, 20)};  RHS = {mp.nstr(rhs, 20)};  |diff| = {mp.nstr(abs(vA-rhs),3)}")
log(f"    -H01 - H0-1 + zeta(2)/2 = {mp.nstr(vB, 20)} (printed in HEXC2 result.md: |diff| = {mp.nstr(abs(vB-rhs),3)})")
assert abs(vA - rhs) < mp.mpf("1e-45") and abs(vB - rhs) > mp.mpf("1e-3")
log("    => certified ID2 has +H_{0,-1}; the display '-H_{0,1}-H_{0,-1}' in HEXC2 result.md is a sign typo")
log("STEP 07 PASS")
with open("progress.log", "a") as f:
    f.write("step_07: my targets == hexc2 targets EXACTLY; cert exact vs my targets; |T_n|<=1e-55; ID2 sign typo found in their result.md\\n")
'''
with open("steps/step_07.py", "w") as f:
    f.write(step07)
import subprocess
r = subprocess.run(["nice", "-n", "10", "<project>/.venv/bin/python", "steps/step_07.py"], capture_output=True, text=True, timeout=3500)
print(r.stdout)
print(r.stderr[-2500:] if r.returncode else "OK")
