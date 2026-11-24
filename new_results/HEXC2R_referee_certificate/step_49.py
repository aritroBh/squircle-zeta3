
step05 = '''"""step_05.py -- referee: exact Q(zeta) arithmetic, Mobius maps, cross-ratios, weight-1 rows.
 (a) my own exact arithmetic in Q(zeta), zeta^2 = zeta - 1, as pairs (a,b) = a + b*zeta;
 (b) enumerate all Mobius maps sending {0,1/2,1,-1,inf} injectively into mu_6 -- exact;
 (c) verify the task map f(z) = (zeta^2 z+1)/(zeta^4 z+1) images -- exact;
 (d) cross-ratio table of the 5-point configuration: which values occur; check
     1/3, -1/3, 2/3, 1/2, -1/2, -2 occur (exact rationals);
 (e) Li4(1/9) = 8(Li4(1/3)+Li4(-1/3)) numeric @200 dps (+ exact series argument);
 (f) weight-1 branch rows vanish numerically at 65 dps (my evaluator);
 (g) family breakdown of the certificate rows."""
import sys, json, time
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "../HEXC2_level6_certificate")
sys.path.insert(0, "../HEXC1_hexagon_certificate")
import mpmath as mp
from fractions import Fraction as Fr
from itertools import product, permutations
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

# ---------- (a) exact Q(zeta): (a,b) = a + b zeta, zeta^2 = zeta - 1 ----------
def qadd(x, y): return (x[0] + y[0], x[1] + y[1])
def qneg(x): return (-x[0], -x[1])
def qsub(x, y): return (x[0] - y[0], x[1] - y[1])
def qmul(x, y):
    a, b = x; c, d = y
    # (a + b z)(c + d z) = ac + (ad+bc) z + bd z^2 = (ac - bd) + (ad + bc + bd) z
    return (a*c - b*d, a*d + b*c + b*d)
def qinv(x):
    a, b = x
    # 1/(a + b z) = (a + b z')/N, z' = 1 - z conjugate; N = a^2 + ab + b^2... compute:
    # (a + b z)(a + b(1-z)) = a^2 + ab z + ab(1-z)... = a^2 + ab + b^2 z(1-z)= a^2+ab+b^2
    n = a*a + a*b + b*b
    assert n != 0
    return ((a + b) / n, -b / n) if False else (Fr(a + b, n), Fr(-b, n))
def qdiv(x, y): return qmul(x, qinv(y))
def qpow(x, n):
    out = (Fr(1), Fr(0))
    for _ in range(n): out = qmul(out, x)
    return out
ZETA = (Fr(0), Fr(1))
MU6 = [qpow(ZETA, k) for k in range(6)]
INF = None
def qmobius(abcd, x):
    """(A z + B)/(C z + D) at x in Q(frac) or INF; returns (num,den) pair or INF."""
    A, B, C, D = abcd
    if x is INF:
        return qdiv(A, C) if C != (0, 0) else INF
    xq = (Fr(x.numerator, x.denominator), Fr(0))
    n = qadd(qmul(A, xq), B); d = qadd(qmul(C, xq), D)
    if d == (0, 0): return INF
    return qdiv(n, d)
# encode rationals into Q(zeta) scalars:
def qr(x):
    x = Fr(x)
    return (x, Fr(0))
# task map: (zeta^2 z + 1)/(zeta^4 z + 1)
TASK = (qpow(ZETA, 2), (Fr(1), Fr(0)), qpow(ZETA, 4), (Fr(1), Fr(0)))
PTS = [Fr(0), Fr(1, 2), Fr(1), Fr(-1), INF]

# ---------- (c) task map images ----------
imgs = [qmobius(TASK, p) for p in PTS]
exp_imgs = []
for im in imgs:
    exp_imgs.append(MU6.index(im) if im in MU6 else None)
log(f"(c) task map images of [0, 1/2, 1, -1, inf]: mu6 exponents = {exp_imgs} (expect [0,1,2,5,4])")
assert exp_imgs == [0, 1, 2, 5, 4]

# ---------- (b) enumerate all Mobius maps sending the 5 points injectively into mu6 ----------
# A Mobius map is determined by the images of 0, 1, inf; it has the form
# f(z) = (a z + b)/(c z + d). Build from images f0, f1, finf in mu6 distinct, then test f(1/2), f(-1).
def mobius_from_3(f0, f1, finf):
    # solve for map with f(0)=f0, f(1)=f1, f(inf)=finf over Q(zeta):
    # f(z) = finf + (f0 - finf) * (f1 - finf) / ( (f1 - finf) - (f0 - finf) z )   when all finite:
    # handle infinities: our candidate images are all in mu6 (finite), so:
    a = qsub(f1, finf); b = qsub(f0, finf)
    # f(z) = (finf * a + b * (f0... )) -- derive: f(z) = (A z + B)/(C z + D) with
    # f(0) = B/D = f0, f(inf) = A/C = finf, f(1) = (A+B)/(C+D) = f1.
    # choose C = 1, A = finf; D free; B = f0 D; (finf + f0 D)/(1 + D) = f1 -> D = (finf - f1)/(f1 - f0)
    C = (Fr(1), Fr(0)); A = finf
    D = qdiv(qsub(finf, f1), qsub(f1, f0))
    B = qmul(f0, D)
    return (A, B, C, D)
found = []
for f0, f1, finf in permutations(MU6, 3):
    M = mobius_from_3(f0, f1, finf)
    ims = [qmobius(M, p) for p in PTS]
    if all((im is not INF and im in MU6) for im in ims) and len(set(ims)) == 5:
        exps = tuple(MU6.index(im) for im in ims)
        found.append(exps)
found = sorted(set(found))
log(f"(b) exactly {len(found)} Mobius maps send {{0,1/2,1,-1,inf}} injectively into mu6 (exact over Q(zeta))")
assert len(found) == 12
theirs = sorted(tuple(m["exps"]) for m in json.load(open("../HEXC2_level6_certificate/mobius_maps.json")))
assert found == theirs, (found, theirs)
log("    my 12 exponent-tuples == mobius_maps.json exactly: " + str(found))
# each map omits exactly one root, 2 maps per root:
from collections import Counter
omit = Counter(tuple(sorted(set(range(6)) - set(f)))[0] for f in found)
log(f"    omitted-root distribution: {dict(omit)} (expect each root omitted twice)")
assert all(v == 2 for v in omit.values())

# ---------- (d) cross-ratio table ----------
def cr(e, p0, p1, pinf):
    """(e,p0;p1,pinf) = (e-p0)(p1-pinf)/((e-pinf)(p1-p0)) over Q, INF handled."""
    def sub(a, b):
        if a is INF and b is INF: return Fr(0)  # not reached for distinct points
        if a is INF: return INF
        if b is INF: return INF
        return a - b
    def div(a, b):
        if a is INF and b is INF: return Fr(1)
        if a is INF: return INF
        if b is INF: return Fr(0)
        return a / b
    f1 = div(sub(e, p0), sub(e, pinf)); f2 = div(sub(p1, pinf), sub(p1, p0))
    if f1 is INF or f2 is INF: return INF
    return f1 * f2
tab = {}
for e, p0, p1, pinf in permutations(PTS, 4):
    z = cr(e, p0, p1, pinf)
    if z is not INF:
        tab.setdefault(z, []).append((e, p0, p1, pinf))
need = [Fr(1, 3), Fr(-1, 3), Fr(2, 3), Fr(1, 2), Fr(-1, 2), Fr(-2)]
log(f"(d) cross-ratio values of the configuration: {sorted(tab.keys())}")
for z in need:
    reps = tab[z]
    log(f"    {z}: {len(reps)} representations, e.g. (e,p0,p1,pinf) = {reps[0]}")
    assert z in tab
assert Fr(1, 9) not in tab
log("    1/9 is NOT a cross-ratio of the configuration (as stated in result.md)")

# ---------- (e) Li4(1/9) distribution ----------
mp.mp.dps = 200
lhs = mp.polylog(4, mp.mpf(1)/9) - 8*(mp.polylog(4, mp.mpf(1)/3) + mp.polylog(4, -mp.mpf(1)/3))
log(f"(e) Li4(1/9) - 8(Li4(1/3)+Li4(-1/3)) = {mp.nstr(lhs, 3)} at 200 dps "
    "(exact: split the absolutely convergent series for Li4(z^2) even/odd)")
assert abs(lhs) < mp.mpf("1e-190")

# ---------- (f) weight-1 branch rows ----------
import reflib as R
mp.mp.dps = 65
v = [R.cmzv_value((j,)) for j in range(1, 6)]
r1 = v[0] + v[4]
r2 = v[0] - v[4] - 2*v[1] + 2*v[3]
log(f"(f) branch rows: (1)+(5) = {mp.nstr(abs(r1),3)}; (1)-(5)-2(2)+2(4) = {mp.nstr(abs(r2),3)}")
assert abs(r1) < mp.mpf("1e-60") and abs(r2) < mp.mpf("1e-60")

# ---------- (g) family breakdown of certificate rows ----------
import cmzv6
cert23 = json.load(open("../HEXC2_level6_certificate/cert_w23.json"))
cert4 = json.load(open("../HEXC2_level6_certificate/cert_w4.json"))
lower = {}
for n in (1, 2, 3, 4):
    rows, counts = cmzv6.make_relations(n, lower)
    lower[n] = rows
    if n == 1: continue
    fam = {}
    def emit(row, f):
        r = cmzv6.primitive(row)
        if r and r not in fam: fam[r] = f
    for a in range(1, n):
        b = n - a
        for u in cmzv6.words(a, False):
            u_adm = u[0] != 0
            if u_adm and a > b: continue
            for vv in cmzv6.words(b):
                if u_adm and a == b and u > vv: continue
                emit(cmzv6.ds(u, vv), 'DS')
    for d in (2, 3, 6):
        for row in cmzv6.distribution(n, d): emit(row, f'dist{d}')
    for k, basis in lower.items():
        if k >= n: continue
        for row in basis:
            for w in cmzv6.words(n - k):
                emit(cmzv6.mul(dict(row), {w: 1}), f'lift{k}')
    prow = cert23[str(n)]["prow"] if n < 4 else cert4["prow"]
    bd = Counter(fam[rows[i]] for i in prow)
    log(f"(g) weight {n} certificate rows by family: {dict(bd)}")
log("STEP 05 PASS")
with open("progress.log", "a") as f:
    f.write("step_05: 12 Mobius maps + task-map images + cross-ratios exact over Q(zeta); "
            "Li4(1/9) dist 1e-190; branch rows 1e-60\\n")
'''
with open("steps/step_05.py", "w") as f:
    f.write(step05)
import subprocess
r = subprocess.run(["nice", "-n", "10", "<project>/.venv/bin/python", "steps/step_05.py"], capture_output=True, text=True, timeout=3000)
print(r.stdout)
print(r.stderr[-2500:] if r.returncode else "OK")
