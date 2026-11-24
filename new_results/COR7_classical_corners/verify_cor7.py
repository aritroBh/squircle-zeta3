"""verify_cor7.py -- COR7: classical corners of M_3(alpha).
Checks (all with printed PASS/FAIL and residuals):
  A. Task 1 (c=1/3): weight split consistency, W2/W3/W4 closed forms vs stored 1145-digit values,
     quadrature vs HPL at dps 60, COR6 T4 Lie class = Li4 part of the stored W4 relation,
     zeta(3) coefficient -26.
  B. Task 2 (census): rational c0 with rank <= 2 are exactly {1/3, 1/2} (enumerated to denom 300);
     the 5 quadratic classical corners: exact letter-span rank 2 and true delta_22 = 0;
     Li4 classes at S-units span Q^3 and contain every component of class(Omega);
     dim-3 controls (sqrt17-4, 1/phi^2, (sqrt2-1)^4, rational 1/5) all have nonzero true delta_22.
  C. Task 3 (c=sqrt2-1): W3/W4 closed forms vs stored, W2 = -Li2(b^2)/2 + pi^2/12,
     Lie class of the Li4 part = exact class(Omega), quadrature vs HPL at dps 60.
  D. COR6 table: recompute the degree <= 2 rows (true vs naive delta_22) and compare with cor6_table.json.
Runs in < 5 min. Requires: cor7_pt.py (workdir), COR6 modules (read-only), stored values (workdir).
"""
import sys, json, time
from fractions import Fraction as Fr

COR7 = "<project>/research/08_theory/COR7_classical_corners"
COR6 = "<project>/research/08_theory/COR6_true_point_symbols"
sys.path.insert(0, COR7)
sys.path.insert(0, COR6)
import mpmath as mp
from cor7_pt import QuadField, QuadPoint
from cor6_sym import W4Quotient, lneg
from cor6_ft import ft_coords, ft_delta22
from cor6_pipe import W4TERMS
from cor6_table import Point

T0 = time.time()
PASS = []
def report(name, ok, extra=""):
    PASS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {extra}", flush=True)

mp.mp.dps = 420
def Li(s, x): return mp.polylog(s, x)

# ============ A. Task 1: c = 1/3 ============
print("=== A. c = 1/3 corner ===")
M3 = mp.mpf(open(COR7 + "/c13_M3_hpl_1150.txt").read())
W4v = mp.mpf(open(COR7 + "/c13_W4_hpl_1150.txt").read())
W3v = mp.mpf(open(COR7 + "/c13_W3_hpl_1150.txt").read())
W2v = mp.mpf(open(COR7 + "/c13_W2_hpl_1150.txt").read())
c = mp.mpf(1)/3
l2, l3, z2, z3, pi2 = mp.log(2), mp.log(3), mp.zeta(2), mp.zeta(3), mp.pi**2
report("A1 W4+W3+W2 = M3 (stored, 1145 digits)", abs(W4v + W3v + W2v - M3) < mp.mpf(10)**(-400),
       mp.nstr(abs(W4v + W3v + W2v - M3), 3))
W2cf = -Li(2, c) - Li(2, -c) + pi2/12
report("A2 W2 = -Li2(1/3) - Li2(-1/3) + pi^2/12", abs(W2v - W2cf) < mp.mpf(10)**(-390),
       mp.nstr(abs(W2v - W2cf), 3))
r3 = json.load(open(COR7 + "/c13_W3_relation.json"))
l313, l31m, l213 = Li(3, c), Li(3, -mp.mpf(1)/2), Li(2, c)
B3 = [z3, l313, l31m, pi2*l2, pi2*l3, l213*l2, l213*l3, l2**3, l2**2*l3, l2*l3**2, l3**3]
W3cf = -sum(mp.mpf(ci)*b for ci, b in zip(r3["coeffs"], B3))/r3["a0"]
report("A3 W3 closed form (identical at dps 600/900 when found)", abs(W3v - W3cf) < mp.mpf(10)**(-390),
       mp.nstr(abs(W3v - W3cf), 3))
report("A3b zeta(3) coefficient is -26", r3["coeffs"][0] == 312 and r3["a0"] == 12)
r4 = json.load(open(COR7 + "/c13_F4_relation.json"))
a0, cf4 = r4["a0"], r4["coeffs"]
B4 = [Li(4, -2), Li(4, c), Li(4, mp.mpf(2)/3), Li(4, mp.mpf(1)/2), Li(4, mp.mpf(1)/9),
      l313*l2, l313*l3, l31m*l2, l31m*l3, z3*l2, z3*l3,
      l213**2, l213*l2**2, l213*l2*l3, l213*l3**2, pi2*l213,
      pi2*l2**2, pi2*l2*l3, pi2*l3**2, mp.pi**4,
      l2**4, l2**3*l3, l2**2*l3**2, l2*l3**3, l3**4]
W4cf = -sum(mp.mpf(ci)*b for ci, b in zip(cf4, B4))/a0
report("A4 W4 closed form (identical at dps 700/1000 when found)", abs(W4v - W4cf) < mp.mpf(10)**(-390),
       mp.nstr(abs(W4v - W4cf), 3))
import cor7_lib as Lq
mp.mp.dps = 60
M3q = Lq.corner_coeffs_hp(2*mp.acos(c), 3)[3]
mp.mp.dps = 420
report("A5 quadrature(dps 60) vs stored HPL M3", abs(M3q - M3) < mp.mpf(10)**(-55),
       mp.nstr(abs(M3q - M3), 3))

# A6: T4 Lie class (exact symbol algebra)
from cor6_sym import sadd, sscale, sym_H
from cor6_table import vecQ
Q2 = W4Quotient(2)
def li4_sym(x, co=Fr(1)):
    vx, vo = vecQ(x), vecQ(1 - x)
    return {(lneg(vo), vx, vx, vx): co}
ST4 = {}
for x, co in ((Fr(-8), Fr(-11, 81)), (Fr(-3), Fr(-2)), (Fr(2, 3), Fr(20, 3))):
    ST4 = sadd(ST4, li4_sym(x, co))
SOm3 = {}
for e, w, val in W4TERMS:
    SOm3 = sadd(SOm3, sscale(val * Fr(1, 3)**e, sym_H(w, Fr(1, 3), vecQ)))
Snum = {}
for x, co in ((Fr(-2), Fr(-22, 3)), (Fr(1, 3), Fr(-17, 2)), (Fr(2, 3), Fr(20, 3)),
              (Fr(1, 2), Fr(110, 9)), (Fr(1, 9), Fr(25, 48))):
    Snum = sadd(Snum, li4_sym(x, co))
ok = Q2.coords(SOm3) == Q2.coords(ST4) == Q2.coords(Snum)
report("A6 class(Omega(1/3)) = COR6 T4 class = Li4 part of numerical W4 (coords (14, 20/3, 8/3))",
       ok, str(Q2.coords(SOm3)))

# ============ B. Task 2: census ============
print("=== B. census ===")
from math import gcd
import sympy as sp
def letter_rank_rational(cc):
    letters = [Fr(2), cc, 1-cc, 1+cc]
    primes = set(); facs = []
    for L in letters:
        f = {}
        for q in (abs(L.numerator), L.denominator):
            for p, e in sp.factorint(q).items():
                f[p] = f.get(p, 0) + (e if q == abs(L.numerator) else -e)
                primes.add(p)
        facs.append(f)
    primes = sorted(primes)
    M = [[f.get(p, 0) for p in primes] for f in facs]
    R = 0
    for cc2 in range(len(M[0])):
        p = next((r for r in range(R, len(M)) if M[r][cc2]), None)
        if p is None: continue
        M[R], M[p] = M[p], M[R]
        q = M[R][cc2]; M[R] = [x/q for x in M[R]]
        for r in range(len(M)):
            if r != R and M[r][cc2]:
                M[r] = [a - b*M[R][cc2] for a, b in zip(M[r], M[R])]
        R += 1
    return R
sols = [(a, b) for b in range(2, 301) for a in range(1, b)
        if gcd(a, b) == 1 and letter_rank_rational(Fr(a, b)) <= 2]
report("B1 rational c0 (denom <= 300) with letter rank <= 2: exactly {1/3, 1/2}",
       set(sols) == {(1, 2), (1, 3)}, str(sols))

CORNERS = [
    (2, "sqrt2-1", (Fr(-1), Fr(1))),
    (2, "(sqrt2-1)^2", (Fr(3), Fr(-2))),
    (2, "1/sqrt2 (square)", (Fr(0), Fr(1, 2))),
    (5, "1/phi (golden)", (Fr(-1, 2), Fr(1, 2))),
    (5, "1/phi^3 = sqrt5-2", (Fr(-2), Fr(1))),
]
QFs = {2: QuadField(2), 5: QuadField(5)}
okB2 = True
for D, nm, c0 in CORNERS:
    P = QuadPoint(QFs[D], nm, c0, primes=[2])
    r, labs, vecs = P.letter_span_rank()
    SOm = P.omega_symbol(W4TERMS)
    d22 = ft_delta22(SOm, 2, len(P.gensE))
    ok = (r == 2) and (len(d22) == 0)
    okB2 = okB2 and ok
    report(f"B2 {nm}: exact letter-span rank 2 and TRUE delta_22 = 0", ok)
def li4_coords(P, x):
    vx, vo = P.vec(x), P.vec(P.F.sub(P.F.one(), x))
    return Q2.coords({(lneg(vo), vx, vx, vx): Fr(1)})
def rankQ(vecs):
    M = [list(v) for v in vecs]; R = 0
    for cc2 in range(len(M[0])):
        p = next((r for r in range(R, len(M)) if M[r][cc2]), None)
        if p is None: continue
        M[R], M[p] = M[p], M[R]
        q = M[R][cc2]; M[R] = [x/q for x in M[R]]
        for r in range(len(M)):
            if r != R and M[r][cc2]:
                M[r] = [a - b*M[R][cc2] for a, b in zip(M[r], M[R])]
        R += 1
    return R
import itertools
okB3 = True
for D, nm, c0 in CORNERS:
    QF = QFs[D]
    P = QuadPoint(QF, nm, c0, primes=[2])
    F = P.F
    gens = [g[0] for g in P.gensE]
    argset = []
    for tup in itertools.product(range(-4, 5), repeat=len(gens)):
        v = F.one()
        for e_, g in zip(tup, gens):
            if e_: v = F.mul(v, F.pow(g, e_))
        for sgn in (1, -1):
            x = v if sgn == 1 else F.neg(v)
            if F.eq(x, F.one()): continue
            if P.vec(F.sub(F.one(), x)) is not None and not any(F.eq(x, y) for y in argset):
                argset.append(x)
    span = [li4_coords(P, x) for x in argset]
    cls = ft_coords(Q2, P.omega_symbol(W4TERMS), 2)
    rspan = rankQ(span)
    inspan = all(rankQ(span + [list(comp)]) == rspan for comp in cls)
    ok = (rspan == 3) and inspan
    okB3 = okB3 and ok
    report(f"B3 {nm}: Li4 classes span Q^3 and contain class(Omega) ({len(argset)} args)", ok)
# dim-3 controls
CTRL = [(17, "sqrt17-4", (Fr(-4), Fr(1)), [2]), (5, "1/phi^2", (Fr(3, 2), Fr(-1, 2)), [2, 5]),
        (2, "(sqrt2-1)^4", (Fr(17), Fr(-12)), [2, 3])]
if 17 not in QFs: QFs[17] = QuadField(17)
for D, nm, c0, pr in CTRL:
    P = QuadPoint(QFs[D], nm, c0, primes=pr)
    r, labs, vecs = P.letter_span_rank()
    d22 = ft_delta22(P.omega_symbol(W4TERMS), 2, len(P.gensE))
    report(f"B4 {nm}: letter-span rank {r} >= 3 and TRUE delta_22 nonzero", r >= 3 and len(d22) > 0)
P15 = Point("c=1/5", [-1, 1], (Fr(1, 5),), basis="auto", verbose=False)
d15 = ft_delta22(P15.omega_symbol(), 1, len(P15.basis))
report("B4 rational c=1/5 (dimV 3): TRUE delta_22 nonzero", len(d15) > 0, f"{len(d15)} comp")

# ============ C. Task 3: c = sqrt2-1 ============
print("=== C. c = sqrt2 - 1 closed form ===")
mp.mp.dps = 420
M3b = mp.mpf(open(COR7 + "/b21_M3_hpl_1150.txt").read())
W4b = mp.mpf(open(COR7 + "/b21_W4_hpl_1150.txt").read())
W3b = mp.mpf(open(COR7 + "/b21_W3_hpl_1150.txt").read())
W2b = mp.mpf(open(COR7 + "/b21_W2_hpl_1150.txt").read())
b = mp.sqrt(2) - 1; rt2 = mp.sqrt(2)
Llog = mp.log(1 + rt2); lam = mp.log(2)
Li3b, Li3mb, Li3s = Li(3, b), Li(3, -b), Li(3, 1/rt2)
Li2b = Li(2, b)
r3b = json.load(open(COR7 + "/b21_W3_relation.json"))
a03, cf3 = r3b["a0"], r3b["coeffs"]
B3b = [z3, Li3b, Li3mb, Li3s, Llog*Li2b, lam*Li2b, pi2*Llog, pi2*lam, lam**3, lam**2*Llog, lam*Llog**2, Llog**3]
nb3 = len(B3b)
W3bcf = -sum(mp.mpf(cf3[i])*B3b[i] + mp.mpf(cf3[i+nb3])*rt2*B3b[i] for i in range(nb3))/a03
report("C1 W3 closed form (identical at dps 700/1000 when found)", abs(W3b - W3bcf) < mp.mpf(10)**(-390),
       mp.nstr(abs(W3b - W3bcf), 3))
r4b = json.load(open(COR7 + "/b21_W4_relation.json"))
a04, cf4b = r4b["a0"], r4b["coeffs"]
B4b = [Li(4, b), Li(4, -b), Li(4, -b*b), Li(4, 1/rt2), Li(4, mp.mpf(1)/2), Li(4, 1-1/rt2), Li(4, 2-rt2),
       Llog*Li3b, Llog*Li3mb, Llog*Li3s, Llog*z3, lam*Li3b, lam*Li3mb, lam*Li3s, lam*z3,
       Li2b**2, Li2b*Llog**2, Li2b*lam*Llog, Li2b*lam**2, pi2*Li2b,
       mp.pi**4, pi2*Llog**2, pi2*lam*Llog, pi2*lam**2,
       Llog**4, Llog**3*lam, Llog**2*lam**2, Llog*lam**3, lam**4]
nb4 = len(B4b)
W4bcf = -sum(mp.mpf(cf4b[i])*B4b[i] + mp.mpf(cf4b[i+nb4])*rt2*B4b[i] for i in range(nb4))/a04
report("C2 W4 closed form (identical at dps 700/1000 when found)", abs(W4b - W4bcf) < mp.mpf(10)**(-390),
       mp.nstr(abs(W4b - W4bcf), 3))
W2bcf = -Li(2, b*b)/2 + pi2/12
report("C3 W2 = -Li2((sqrt2-1)^2)/2 + pi^2/12 (duplication)", abs(W2b - W2bcf) < mp.mpf(10)**(-390),
       mp.nstr(abs(W2b - W2bcf), 3))
report("C4 M3 = W4+W3+W2 closed forms", abs(M3b - (W4bcf + W3bcf + W2bcf)) < mp.mpf(10)**(-380),
       mp.nstr(abs(M3b - (W4bcf + W3bcf + W2bcf)), 3))
# Lie class match (exact)
P = QuadPoint(QFs[2], "b", (Fr(-1), Fr(1)), primes=[2])
argsC = {"b": (Fr(-1), Fr(1)), "-b": (Fr(1), Fr(-1)), "-b^2": (Fr(-3), Fr(2)),
         "1/sqrt2": (Fr(0), Fr(1, 2)), "1/2": (Fr(1, 2), Fr(0)),
         "1-1/sqrt2": (Fr(1), Fr(-1, 2)), "2-sqrt2": (Fr(2), Fr(-1))}
coefC = {"b": (Fr(22), Fr(-31, 2)), "-b": (Fr(24), Fr(-27, 2)), "-b^2": (Fr(-5, 2), Fr(13, 8)),
         "1/sqrt2": (Fr(0), Fr(-1)), "1/2": (Fr(3, 4), Fr(0)),
         "1-1/sqrt2": (Fr(-4), Fr(4)), "2-sqrt2": (Fr(-2), Fr(3))}
comp0 = [Fr(0)]*3; comp1 = [Fr(0)]*3
for nm, x in argsC.items():
    v = li4_coords(P, x)
    p, q = coefC[nm]
    for j in range(3):
        comp0[j] += p*v[j]; comp1[j] += q*v[j]
cls = ft_coords(Q2, P.omega_symbol(W4TERMS), 2)
report("C5 Li4 part of W4 closed form = exact class(Omega(sqrt2-1)) ((8,-2,-6),(-3,1,8))",
       tuple(comp0) == tuple(cls[0]) and tuple(comp1) == tuple(cls[1]))
mp.mp.dps = 60
M3q2 = Lq.corner_coeffs_hp(2*mp.acos(mp.sqrt(2)-1), 3)[3]
mp.mp.dps = 420
report("C6 quadrature(dps 60) vs stored HPL M3", abs(M3q2 - M3b) < mp.mpf(10)**(-55),
       mp.nstr(abs(M3q2 - M3b), 3))

# ============ D. COR6 table (degree <= 2 rows recomputed) ============
print("=== D. COR6 table (deg <= 2 rows) ===")
tab = json.load(open("<project>/research/08_theory/COR6_true_point_symbols/cor6_table.json"))
expect = {row["label"]: row for row in tab}
ROWS = [
    ("n=3 triangle c=sqrt3/2", [-3, 0, 1], (Fr(0), Fr(1, 2))),
    ("n=4 square c=1/sqrt2", [-2, 0, 1], (Fr(0), Fr(1, 2))),
    ("n=6 hexagon c=1/2", [-1, 1], (Fr(1, 2),)),
    ("n=10 c=sin(pi/10)", [-5, 0, 1], (Fr(-1, 4), Fr(1, 4))),
    ("c=1/3", [-1, 1], (Fr(1, 3),)),
    ("c=2/5", [-1, 1], (Fr(2, 5),)),
    ("c=3/7", [-1, 1], (Fr(3, 7),)),
]
for label, mpoly, c0 in ROWS:
    Pt = Point(label, mpoly, c0, basis="auto", verbose=False)
    dT = ft_delta22(Pt.omega_symbol(), Pt.F.d, len(Pt.basis))
    dN = Pt.naive_delta22()
    row = expect[label]
    ok = (len(dT) == 0) == row["true_zero"] and (len(dN) == 0) == row["naive_zero"] and (dT != dN) == row["differ"]
    report(f"D {label}: true_zero={row['true_zero']} naive_zero={row['naive_zero']} differ={row['differ']} reproduced", ok)

print(f"\nTOTAL: {sum(1 for _, ok in PASS if ok)}/{len(PASS)} PASS, time {time.time()-T0:.0f}s")
if not all(ok for _, ok in PASS):
    sys.exit(1)
