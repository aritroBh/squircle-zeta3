verify = r'''"""verify_cor13c.py -- COR13C certificate verification. Runtime target < 3 min.

PART A (the proof, exact rational arithmetic):
  A1. Exact re-enumeration: for N in {6,12,18,24,30,36} the number of injective Mobius maps
      {0,1/3,1,-1,inf} -> mu_N over Q(zeta_N) is 2N (mod-p filter with sharp primes
      [(p-1)/N >= 50; complete superset] + exact verification in Q(zeta_N)); at N = 6 also
      a direct exact enumeration over Q(zeta_6) (HEXC2 step_18 style). No map needs 0 or inf.
  A2. For n in {2,3,4}: rebuild the proved level-6 relation rows with cmzv6.make_relations
      (regularized double shuffle [Ihara-Kaneko-Zagier], Zhao regularized distribution
      d in {2,3,6}, weight-1 branch identities, shuffle lifts), rebuild the target vector
      T_n (ID_n of COR7's c0=1/3 closed form: corner side minus closed form, converted to
      level-6 CMZV words via the Mobius map), load the certificate (row indices prow,
      rational multipliers y) from cert13_w23.json / cert13_w4.json, and check
          sum_i y_i rows[prow_i] - T_n == 0   exactly, on every admissible word.
      Then T_n is an exact rational combination of proved relations, hence evaluates to 0.
  A3. Bonus: zeta(2) = pi^2/6 itself is an exact combination of the weight-2 proved rows.

PART B (numeric audit of the conversion code, 62 dps; not part of the proof):
  the 12 Mobius map images, all 38 HPL-word conversions at 1/3 vs H_series, all 9 Li atoms
  vs polylog, the 7 elementary atoms, the three target values, the Li4(1/9) distribution
  input, and a cross-check of the two independent CMZV evaluators.
"""
import os, sys, json, time, re
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
THEORY = "<project>/research/08_theory"
sys.path.insert(0, THEORY + "/HEXC1_hexagon_certificate")
sys.path.insert(0, THEORY + "/COR7_classical_corners")
from fractions import Fraction as Fr
T0 = time.monotonic()
def log(s):
    print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

PART_B = ("--fast" not in sys.argv)

# ---------------- Part A ----------------
log("PART A: exact certificate verification")

# A1: exact map enumeration
from itertools import permutations
import sympy as sp
from cyclo import Cyclo, mobius_images
def prime_for(N, kmin=50):
    k = kmin
    while True:
        p = k*N + 1
        if sp.isprime(p): return int(p)
        k += 1
okA = True
counts = {}
for N in (6, 12, 18, 24, 30, 36):
    p = prime_for(N)
    g = int(sp.primitive_root(p)); z = pow(g, (p-1)//N, p)
    roots = [pow(z, k, p) for k in range(N)]; rset = set(roots)
    cand = []
    for i, j, k in permutations(range(N), 3):
        a, b, c = roots[i], roots[j], roots[k]
        bma = (b-a) % p; cmb = (c-b) % p
        t1 = c*bma % p; t2 = a*cmb % p
        cls = []
        for num, den in (((t1 + 3*t2) % p, (bma + 3*cmb) % p), ((t2 - t1) % p, (cmb - bma) % p)):
            if den == 0: cls.append('inf')
            elif num == 0: cls.append('zero')
            else:
                cls.append('root' if num*pow(den, -1, p) % p in rset else 'other')
        if 'other' not in cls:
            cand.append((i, j, k))
    C = Cyclo(N)
    strict = 0; extra = 0
    for i, j, k in cand:
        im = list(mobius_images(C, i, j, k))
        if all(x[0] == 'root' for x in im):
            exps = (i, im[0][1], j, im[1][1], k)
            if len(set(exps)) == 5: strict += 1
        elif all(x[0] in ('root', 'zero', 'inf') for x in im):
            extra += 1
    counts[N] = strict
    good = (strict == 2*N and extra == 0)
    okA = okA and good
    log(f"  [{'PASS' if good else 'FAIL'}] N={N}: {strict} injective maps into mu_N "
        f"(expected 2N={2*N}), maps using 0/inf: {extra}")
# direct exact enumeration at N=6 (HEXC2 step_18 style, Fractions over Q(zeta_6))
from fractions import Fraction as Fq
def zadd(x,y): return (x[0]+y[0], x[1]+y[1])
def zneg(x): return (-x[0], -x[1])
def zsub(x,y): return (x[0]-y[0], x[1]-y[1])
def zmul(x,y):
    a,b=x; c,d=y
    return (a*c-b*d, a*d+b*c+b*d)
def zinv(x):
    a,b=x; n = a*a+a*b+b*b
    return ((a+b)/n, -b/n)
def zdiv(x,y): return zmul(x, zinv(y))
Z0=(Fq(0),Fq(0)); Z1=(Fq(1),Fq(0)); ZETA=(Fq(0),Fq(1))
MU6=[]
r=Z1
for _ in range(6): MU6.append(r); r=zmul(r,ZETA)
INF=None
def cr_mobius(a,b,c):
    if a is INF: return [Z0, zsub(b,c), Z1, zneg(c)]
    if c is INF: return [Z1, zneg(a), Z0, zsub(b,a)]
    bc = zsub(b,c); ba = zsub(b,a)
    return [bc, zneg(zmul(a,bc)), ba, zneg(zmul(c,ba))]
def minv(M):
    A,B,C,D = M; return [D, zneg(B), zneg(C), A]
def mcomp(M,Nn):
    A,B,C,D = M; a,b,c,d = Nn
    return [zadd(zmul(A,a),zmul(B,c)), zadd(zmul(A,b),zmul(B,d)),
            zadd(zmul(C,a),zmul(D,c)), zadd(zmul(C,b),zmul(D,d))]
def mapply(M, z):
    A,B,C,D = M
    if z is INF:
        return INF if C==Z0 else zdiv(A,C)
    num = zadd(zmul(A,z),B); den = zadd(zmul(C,z),D)
    return INF if den==Z0 else zdiv(num,den)
n6 = 0; exps6 = []
for w0,w1,wI in permutations(range(6),3):
    M = mcomp(minv(cr_mobius(MU6[w0],MU6[w1],MU6[wI])), cr_mobius(Z0,Z1,INF))
    imgs = [mapply(M,p) for p in [Z0,(Fq(1,3),Fq(0)),Z1,(Fq(-1),Fq(0)),INF]]
    if all(im is not INF and im in MU6 for im in imgs):
        exps = tuple(MU6.index(im) for im in imgs)
        if len(set(exps))==5:
            n6 += 1; exps6.append(exps)
stored = [tuple(m['exps']) for m in json.load(open(HERE + "/mobius_maps_13.json"))]
good = (n6 == 12 and sorted(exps6) == sorted(stored))
okA = okA and good
log(f"  [{'PASS' if good else 'FAIL'}] N=6 direct exact enumeration: {n6} maps, "
    f"matches mobius_maps_13.json: {sorted(exps6)==sorted(stored)}")

# A2: certificates
import cmzv6
from cor13c_common import target2, target3, target4, cross_table
from flint import fmpq_mat
# cross-ratio availability of all Li arguments used
tab = cross_table()
need = [Fq(1,3), Fq(-1,3), Fq(2,3), Fq(1,2), Fq(-1,2), Fq(-2)]
good = all(zz in tab for zz in need)
okA = okA and good
log(f"  [{'PASS' if good else 'FAIL'}] all 6 Li arguments are cross-ratios of the configuration")

cert23 = json.load(open(HERE + "/cert13_w23.json"))
cert4 = json.load(open(HERE + "/cert13_w4.json"))
def exact_verify(n, rows, cols, prow, y, T):
    sel = [rows[i] for i in prow]
    ids = {w: i for i, w in enumerate(cols)}
    assert set(T) <= set(cols), "target word outside admissible words"
    S = fmpq_mat(len(sel), len(cols))
    for i, row in enumerate(sel):
        for w, c in row:
            S[i, ids[w]] = c
    Y = fmpq_mat(len(sel), 1)
    for i, c in enumerate(y):
        Y[i, 0] = fmpq_mat(1, 1, [str(c)])[0, 0]
    res = Y.transpose() * S
    bad = 0
    for w, j in ids.items():
        if Fr(str(res[0, j])) != T.get(w, Fr(0)):
            bad += 1
    return bad
lower = {}
for n in (1, 2, 3, 4):
    rows, counts = cmzv6.make_relations(n, lower)
    lower[n] = rows
    log(f"  built relations weight {n}: {len(rows)} rows, families {counts}")
    if n < 2: continue
    T = {2: target2, 3: target3, 4: target4}[n]()
    cols = cmzv6.words(n)
    if n < 4:
        C = cert23[str(n)]
        prow, y = C["prow"], [Fr(int(a), int(b)) for a, b in C["y"]]
    else:
        prow, y = cert4["prow"], [Fr(s) for s in cert4["y"]]
    assert [list(w) for w in cols] == (C["cols"] if n < 4 else cert4["cols"])
    bad = exact_verify(n, rows, cols, prow, y, T)
    okA = okA and (bad == 0)
    log(f"  [{'PASS' if bad==0 else 'FAIL'}] weight {n}: certificate rows={len(prow)} "
        f"target_support={len(T)} residual_words={bad}")
    if n == 4:
        c1111 = T.get((1,1,1,1), Fr(0))
        log(f"        (pi^4 comes out exactly: T4 coefficient of (1,1,1,1) = {c1111}, killed by certificate)")
# A3: Euler
euler = json.load(open(HERE + "/cert13_euler.json"))
ET = {(-1, 0): Fr(1), (1, 1): Fr(3)}
cols2 = cmzv6.words(2)
assert [list(w) for w in cols2] == euler["cols"]
bad = exact_verify(2, lower[2], cols2, euler["prow"], [Fr(int(a), int(b)) for a, b in euler["y"]], ET)
okA = okA and (bad == 0)
log(f"  [{'PASS' if bad==0 else 'FAIL'}] bonus: zeta(2)=pi^2/6 exact ({len(euler['prow'])} rows), residual_words={bad}")
log("PART A " + ("PASSED" if okA else "FAILED"))

if not PART_B:
    log("(--fast: skipping numeric audit)")
    raise SystemExit(0 if okA else 1)

# ---------------- Part B ----------------
log("PART B: numeric audit at 62 dps (validation of the conversion code)")
import mpmath as mp
import fast_cmzv as FC
import numeric_cmzv as NC
from cor13c_common import (hpl_to_cmzv, li_atom, atom_log2, atom_log3, atom_ipi,
                           atom_zeta2, atom_zeta3, atom_pi2, atom_pi4, MAPS, _fit_terms)
from cor7_lib import H_series
mp.mp.dps = 65
FC.setup(62); FC.finish_setup()
def combo(d):
    tot = mp.mpc(0)
    for w, c in d.items():
        v = FC.cmzv_value(w)
        tot += mp.mpf(c.numerator) / c.denominator * mp.mpc(v.real.mid(), v.imag.mid())
    return tot
okB = True
def check(name, resid, tol):
    global okB
    good = abs(resid) < tol
    okB = okB and good
    log(f"  [{'PASS' if good else 'FAIL'}] {name}: |resid| = {mp.nstr(abs(resid), 3)}")
# B0 maps
zeta = mp.e ** (mp.pi * mp.j / 3)
maperr = 0
for exps in MAPS:
    e0, eth, e1, em, einf = exps
    a, b, c = zeta**e0, zeta**e1, zeta**einf
    for zval, e in [(mp.mpf(1)/3, eth), (mp.mpf(-1), em)]:
        got = (c*(b-a)*zval + a*(c-b)) / ((b-a)*zval + (c-b))
        maperr = max(maperr, abs(got - zeta**e))
check("12 Mobius maps on {1/3,-1} (numeric)", maperr, mp.mpf("1e-40"))
# B1 HPL words
WORDS = sorted({w for name in ["hpl_fit_k3_(0,0,0).json", "hpl_fit_k3_(1,0,0).json"]
                for e, w, v in _fit_terms(name) if w}, key=lambda x: (len(x), x))
worst = 0; wword = None
for w in WORDS:
    err = abs(combo(hpl_to_cmzv(w)) - H_series(w, mp.mpf(1)/3))
    if err > worst: worst, wword = err, w
check(f"all {len(WORDS)} HPL-word conversions at 1/3 (worst: H{list(wword)})", worst, mp.mpf("1e-58"))
# B2 Li atoms
LI = [(2, "1/3"), (2, "-1/3"), (3, "1/3"), (3, "-1/2"), (4, "1/3"), (4, "-1/3"),
      (4, "2/3"), (4, "1/2"), (4, "-2")]
worst = 0; latom = None
for n, z in LI:
    err = abs(combo(li_atom(n, Fr(z))) - mp.polylog(n, mp.mpf(z)))
    if err > worst: worst, latom = err, (n, z)
check(f"all {len(LI)} Li-atom conversions (worst: Li_{latom[0]}({latom[1]}))", worst, mp.mpf("1e-58"))
# B3 elementary atoms
check("zeta(2) atom", combo(atom_zeta2()) - mp.zeta(2), mp.mpf("1e-58"))
check("zeta(3) atom", combo(atom_zeta3()) - mp.zeta(3), mp.mpf("1e-58"))
check("log2 atom", combo(atom_log2()) - mp.log(2), mp.mpf("1e-58"))
check("log3 atom", combo(atom_log3()) - mp.log(3), mp.mpf("1e-58"))
check("i*pi atom", combo(atom_ipi()) - mp.pi * mp.j, mp.mpf("1e-58"))
check("(i*pi)^2 atom = -pi^2", combo(atom_pi2()) + mp.pi**2, mp.mpf("1e-58"))
check("(i*pi)^4 atom = pi^4", combo(atom_pi4()) - mp.pi**4, mp.mpf("1e-58"))
# B4 targets
for name, T in [("ID2", target2), ("ID3", target3), ("ID4", target4)]:
    v = combo(T())
    check(f"target {name} value (support {len(T())})", v, mp.mpf("1e-55"))
# B5 distribution input
lhs = mp.polylog(4, mp.mpf(1)/9) - 8*(mp.polylog(4, mp.mpf(1)/3) + mp.polylog(4, -mp.mpf(1)/3))
check("Li4(1/9) = 8(Li4(1/3)+Li4(-1/3)) (input theorem)", lhs, mp.mpf("1e-58"))
# B6 two evaluators
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
check("two independent CMZV evaluators (12 random words)", worst, mp.mpf("1e-55"))
log("PART B " + ("PASSED" if okB else "FAILED"))
log("VERDICT: " + ("COR13C CERTIFICATE VALID - ID2, ID3, ID4 PROVED" if (okA and okB) else "FAILURE"))
assert okA and okB
'''
open('verify_cor13c.py','w').write(verify)
print("written verify_cor13c.py", len(verify), "bytes")
import subprocess, time
t0 = time.time()
r = subprocess.run(['<project>/.venv/bin/python', 'verify_cor13c.py'], capture_output=True, text=True)
print(r.stdout); print(r.stderr[-3000:] if r.stderr else "")
print(f"elapsed {time.time()-t0:.1f}s")
