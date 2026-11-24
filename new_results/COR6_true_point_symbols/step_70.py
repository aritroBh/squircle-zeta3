verify = r'''"""verify_cor6.py -- COR6: the TRUE symbol of M_3(alpha) at a given angle.
Controls (a)-(d) plus the full corner table (naive vs TRUE delta_22).
Runs < 5 min. Requires: cor6_sym.py, cor6_nf.py, cor6_sunit.py, cor6_ft.py, cor6_table.py.
"""
import sys, json, time
from fractions import Fraction as Fr
sys.path.insert(0, ".")
from cor6_nf import NF
from cor6_sunit import exponent_vector, LV
from cor6_sym import *
from cor6_ft import *
from cor6_table import Point, W4TERMS, vecQ

BASE = "<project>/research/"
t_start = time.time()
PASS = []
def report(name, ok, extra=""):
    PASS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {extra}")

# ---------------- control (a): single words at x = 1/2 ----------------
Q2 = W4Quotient(2)
S = sym_H((0, 0, -1, 1), Fr(1, 2), vecQ)
ok1 = Q2.coords(S) == (Fr(-6), Fr(0), Fr(0))
report("(a1) class H_{0,0,-1,1}(1/2) = 6*class(Li4(-2))", ok1, str(Q2.coords(S)))
ok2 = Q2.coords(sym_H((0, 0, 1, 1), Fr(1, 2), vecQ)) == (Fr(0),)*3
ok3 = Q2.coords(sym_H((0, 1, 0, 1), Fr(1, 2), vecQ)) == (Fr(0),)*3
report("(a2) class Li_{3,1}(1/2,1) = 0", ok2)
report("(a3) class Li_{2,2}(1/2,1) = 0", ok3)
ok4 = True
for xs in ["1/2", "1/3", "2/3", "-2", "1/9", "3/4"]:
    x = Fr(xs)
    Sx = sym_H((0, 0, 0, 1), x, vecQ)
    vx, vo = vecQ(x), vecQ(1 - x)
    SLi = {(lneg(vo), vx, vx, vx): Fr(1)}
    if expand_sym(Sx, 2) != expand_sym(SLi, 2): ok4 = False
report("(a4) S(H_{0,0,0,1}(x)) = S(Li4(x)) exact tensors, 6 x", ok4)
# COR5 word_classes.pkl ground truth (24 words)
import pickle
wc = pickle.load(open(BASE + "08_theory/COR5_corner_classicality/word_classes.pkl", "rb"))
ok5 = all(Q2.coords(sym_H(w, Fr(1, 2), vecQ)) == wc[w][0] for w in wc)
report("(a5) all 24 numeric word classes of COR5 reproduced", ok5)

# ---------------- control (b): whole Omega(1/2) ----------------
SOm = {}
for e, w, val in W4TERMS:
    SOm = sadd(SOm, sscale(val * Fr(1, 2)**e, sym_H(w, Fr(1, 2), vecQ)))
SE = expand_sym(SOm, 2)
T = {}
for w, c in SE.items():
    iw = tuple(0 if l == (Fr(1), Fr(0)) else 1 for l in w)
    T[iw] = T.get(iw, 0) + c
rOm = rho(T)
reading = (rOm.get((0, 1, 0, 0), 0), rOm.get((0, 1, 0, 1), 0), rOm.get((0, 1, 1, 0), 0))
ok6 = reading == (Fr(-55, 8), Fr(13, 8), Fr(13, 8)) and delta22(SOm, 2) == {}
report("(b1) rho(Omega(1/2)) = (-55/8,13/8,13/8) (planner reading)", ok6, str(reading))
def li4_sym(x, co=Fr(1)):
    vx, vo = vecQ(x), vecQ(1 - x)
    return {(lneg(vo), vx, vx, vx): co}
def rho_of(S):
    T = {}
    for w, c in expand_sym(S, 2).items():
        iw = tuple(0 if l == (Fr(1), Fr(0)) else 1 for l in w)
        T[iw] = T.get(iw, 0) + c
    return rho(T)
SCF = {}
for x, co in ((Fr(-2), Fr(-21)), (Fr(1, 3), Fr(-65, 8)), (Fr(2, 3), Fr(13, 2)),
              (Fr(1, 2), Fr(137, 6)), (Fr(1, 9), Fr(39, 64))):
    SCF = sadd(SCF, li4_sym(x, co))
ok7 = rho_of(SCF) == rOm
report("(b2) class(Omega(1/2)) = class(Li4 part of COR5 closed form) (full rho dict)", ok7)
S3 = {}
for x, co in ((Fr(-2), Fr(-21)), (Fr(1, 3), Fr(13, 2)), (Fr(2, 3), Fr(13, 2))):
    S3 = sadd(S3, li4_sym(x, co))
ok8 = rho_of(S3) == rOm
report("(b3) class(Omega(1/2)) = class(-21 Li4(-2) + 13/2 Li4(1/3) + 13/2 Li4(2/3))", ok8)

# ---------------- control (c): square corner ----------------
sq2f = NF([-2, 0, 1])
C4 = Point("c=1/sqrt2 (square, n=4)", [-2, 0, 1], (Fr(0), Fr(1, 2)))
SOm4 = C4.omega_symbol()
d22_4 = ft_delta22(SOm4, 2, 2)
one4 = C4.one
b4 = C4.F.sub(C4.F.theta(), one4)
RLi4 = {}
def li4f(q, x):
    vx, vo = C4.vec(x), C4.vec(C4.F.sub(one4, x))
    return {(lneg(vo), vx, vx, vx): q}
for q, x in [((Fr(1), Fr(-573)), b4), ((Fr(-1), Fr(-573)), C4.F.neg(b4)),
             ((Fr(0), Fr(56)), C4.F.neg(C4.F.mul(b4, b4))),
             ((Fr(0), Fr(17)), C4.c0), ((Fr(-37, 4), Fr(-17, 16)), C4.F.scale(Fr(1, 2), one4)),
             ((Fr(10), Fr(86)), C4.F.sub(one4, C4.c0)),
             ((Fr(-10), Fr(86)), C4.F.sub(C4.F.from_int(2), C4.F.theta()))]:
    RLi4 = ft_add(RLi4, li4f(q, x))
Q4sq = W4Quotient(2)
ok9 = d22_4 == {}
ok10 = ft_coords(Q4sq, SOm4, 2) == ft_coords(Q4sq, RLi4, 2)
report("(c1) delta_22(Omega(1/sqrt2)) = 0", ok9)
report("(c2) class(Omega(1/sqrt2)) = class(Li4 part of thm:W4R)", ok10)

# ---------------- control (d): shuffles ----------------
import random
random.seed(7)
okd = True
for t in range(40):
    lu = random.choice([1, 2]); lv_ = 3 - lu
    u = tuple(random.choice([-1, 0, 1]) for _ in range(lu))
    v = tuple(random.choice([-1, 0, 1]) for _ in range(lv_))
    lhs = sshuffle(sym_H(u, Fr(1, 2), vecQ), sym_H(v, Fr(1, 2), vecQ))
    rhs = {}
    for w in shuffles(u, v):
        rhs = sadd(rhs, sym_H(w, Fr(1, 2), vecQ))
    if expand_sym(lhs, 2) != expand_sym(rhs, 2): okd = False
report("(d) 40 random shuffle controls at c=1/2 give symbol equality", okd)

# ---------------- the table ----------------
print("\n=== TABLE ===")
POINTS = [
    ("n=3 triangle c=sqrt3/2", [-3, 0, 1], (Fr(0), Fr(1, 2)),
     [(Fr(1), Fr(1)), (Fr(0), Fr(1)), (Fr(2), Fr(1))], "auto"),
    ("n=4 square c=1/sqrt2", [-2, 0, 1], (Fr(0), Fr(1, 2)), None, "auto"),
    ("n=5 c=sin(pi/5)", [5, 0, -5, 0, 1], (Fr(0), Fr(1, 2)), None, "auto"),
    ("n=6 hexagon c=1/2", [-1, 1], (Fr(1, 2),), None, "auto"),
    ("n=7 c=sin(pi/7)", [-7, 0, 14, 0, -7, 0, 1], (Fr(0), Fr(1, 2)) + (Fr(0),) * 4, None, "letter-span"),
    ("n=8 octagon c=sin(pi/8)=cos(3pi/8)", [2, 0, -4, 0, 1], (Fr(0), Fr(1, 2)), None, "auto"),
    ("n=9 c=sin(pi/9)", [-3, 0, 9, 0, -6, 0, 1], (Fr(0), Fr(1, 2)) + (Fr(0),) * 4, None, "letter-span"),
    ("n=10 c=sin(pi/10)", [-5, 0, 1], (Fr(-1, 4), Fr(1, 4)), None, "auto"),
    ("n=11 c=sin(pi/11)", [-11, 0, 55, 0, -77, 0, 44, 0, -11, 0, 1], (Fr(0), Fr(1, 2)) + (Fr(0),) * 8, None, "letter-span"),
    ("n=12 c=sin(pi/12)", [1, 0, -4, 0, 1], (Fr(0), Fr(1, 2)), None, "auto"),
    ("c=cos(pi/12) (alpha=pi/6)", [1, 0, -4, 0, 1], None, None, "auto"),
    ("c=1/3", [-1, 1], (Fr(1, 3),), None, "auto"),
    ("c=2/5", [-1, 1], (Fr(2, 5),), None, "auto"),
    ("c=3/7", [-1, 1], (Fr(3, 7),), None, "auto"),
]
table = []
SHARED = {}
for label, mpoly, c0, basis, mode in POINTS:
    t0 = time.time()
    if c0 is None:  # cos(pi/12) = 1/(2 theta) in n=12 field
        c0 = SHARED["c12c"]
    if mode == "auto" and mpoly == [1, 0, -4, 0, 1] and "C12" in SHARED:
        P = Point(label, mpoly, c0, basis=SHARED["C12"].basis)
    else:
        P = Point(label, mpoly, c0, basis=("letter-span" if mode == "letter-span" else (basis if basis else "auto")))
    if label.startswith("n=12"):
        SHARED["C12"] = P
        SHARED["c12c"] = P.F.scale(Fr(1, 2), P.F.inv(P.F.theta()))
    SOmP = P.omega_symbol()
    dT = ft_delta22(SOmP, P.F.d, len(P.basis))
    dN = P.naive_delta22()
    ldim = P.letter_span_dim()
    row = {"label": label, "deg": P.F.d, "dimV": P.dimV, "letter_span_dim": ldim,
           "true_d22": {str(k): [str(x) for x in v] for k, v in sorted(dT.items())},
           "naive_d22": {str(k): [str(x) for x in v] for k, v in sorted(dN.items())},
           "true_zero": len(dT) == 0, "naive_zero": len(dN) == 0,
           "differ": dT != dN, "time": round(time.time() - t0, 1)}
    table.append(row)
    print(f"{label}: deg {P.F.d}, dimV {P.dimV}, letter-span {ldim}, "
          f"TRUE d22 {'ZERO' if not dT else f'{len(dT)} comp'}, naive d22 {'ZERO' if not dN else f'{len(dN)} comp'}, "
          f"differ={dT != dN} ({row['time']}s)")
    # shuffle controls in this field (10 each)
    oks = True
    for t in range(10):
        lu = random.choice([1, 2]); lv_ = 3 - lu
        u = tuple(random.choice([-1, 0, 1]) for _ in range(lu))
        v = tuple(random.choice([-1, 0, 1]) for _ in range(lv_))
        lhs = sshuffle(P.sym_H(u), P.sym_H(v))
        rhs = {}
        for w in shuffles(u, v):
            rhs = sadd(rhs, P.sym_H(w))
        dim = len(P.basis)
        if expand_sym(lhs, dim) != expand_sym(rhs, dim):
            oks = False
    report(f"(d) 10 shuffle controls at {label}", oks)

json.dump(table, open("cor6_table.json", "w"), indent=1)

# c=1/3 classicality
SOm3 = {}
for e, w, val in W4TERMS:
    SOm3 = sadd(SOm3, sscale(val * Fr(1, 3)**e, sym_H(w, Fr(1, 3), vecQ)))
cls3 = Q2.coords(SOm3)
Sc = {}
for x, co in ((Fr(-8), Fr(-11, 81)), (Fr(-3), Fr(-2)), (Fr(2, 3), Fr(20, 3))):
    Sc = sadd(Sc, li4_sym(x, co))
ok11 = Q2.coords(Sc) == cls3
report("(T4) c=1/3: class(Omega) = class(-11/81 Li4(-8) - 2 Li4(-3) + 20/3 Li4(2/3))", ok11, str(cls3))

print(f"\nTOTAL: {sum(1 for _, ok in PASS if ok)}/{len(PASS)} PASS, time {time.time()-t_start:.0f}s")
if not all(ok for _, ok in PASS):
    sys.exit(1)
'''
with open("verify_cor6.py", "w") as f:
    f.write(verify)
import subprocess, sys
r = subprocess.run([sys.executable, "verify_cor6.py"], capture_output=True, text=True, timeout=300)
print(r.stdout)
print("STDERR:", r.stderr[-2000:] if r.stderr else "")
