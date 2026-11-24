"""verify_cor9.py -- COR9 verification script (< 5 min).
Checks: A) machinery (class numbers, fund units, norm-l elements);
B) Task 1 (c0 = (sqrt2-1)^2): HPL values vs stored (dps 400), quadrature dps 60, W2/W3/W4 closed
   forms from stored JSON relations, Lie class = COR7's printed class, assembly;
C) Task 2 (golden corners): same for gol1/gol3, plus the exact weight-3/4 basis relations;
D) Task 3: true_rank on the 5 rank-2 points + controls, delta_22 = 0 (2 points spot), the cor7_pt
   letter_span_rank bug example, census filter sanity on D = 2.
"""
import sys, json, time, re
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms'
sys.path.insert(0, WD)
import mpmath as mp
from fractions import Fraction as Fr

T0 = time.time()
PASS = []
def report(name, ok, extra=""):
    PASS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {extra}", flush=True)

# ============ A. machinery ============
print("=== A. machinery ===")
from cor9_lib import class_number, fund_unit, norm_l_element, L7, QuadField, QuadPoint, li4_coords, omega_class, ft_delta22, W4TERMS
report("A1 class numbers: h(2)=h(46)=h(103)=1, h(10)=2, h(79)=3",
       class_number(2) == 1 and class_number(46) == 1 and class_number(103) == 1
       and class_number(10) == 2 and class_number(79) == 3)
report("A2 norm-2 element in Q(sqrt103): 477^2 - 103*47^2 = 2", 477*477 - 103*47*47 == 2)
e2 = fund_unit(2); e5 = fund_unit(5)
report("A3 fund units eps2 = 1+s2, eps5 = (1+s5)/2",
       e2 == (Fr(1), Fr(1)) and e5 == (Fr(1,2), Fr(1,2)))
el, stat = norm_l_element(17, 2)
report("A4 norm-2 element in Q(sqrt17) exists (2 splits, h=1)", el is not None and abs(el[0]**2 - 17*el[1]**2) == 2)

mp.mp.dps = 420
Li = lambda s, x: mp.polylog(s, x)
z2, z3, pi2 = mp.zeta(2), mp.zeta(3), mp.pi**2

def basis_val(nm, phi, r5, L, lam):
    v = mp.mpf(1)
    for fac in nm.split("*"):
        m = re.match(r"^(Li[234])\((.*)\)(?:\^(\d+))?$", fac)
        if m:
            e = m.group(2).replace("^", "**").replace("sqrt5", "r5")
            e = re.sub(r"(\d)\s*phi", r"\1*phi", e)
            e = re.sub(r"r5\s+phi", "r5*phi", e)
            s_ = int(m.group(1)[2]); x = eval(e, {"r5": r5, "phi": phi})
            vv = Li(s_, x)
            if m.group(3): vv = vv**int(m.group(3))
        else:
            m2 = re.match(r"^(zeta3|pi2|pi4|L|lam)(?:\^(\d+))?$", fac)
            vv = {"zeta3": z3, "pi2": pi2, "pi4": mp.pi**4, "L": L, "lam": lam}[m2.group(1)]
            if m2.group(2): vv = vv**int(m2.group(2))
        v *= vv
    return v

def check_relation(tag, wf, basis_vals_list, rt):
    r = json.load(open(WD + f"/{tag}_{wf}_relation.json"))
    a0, cf = r["a0"], r["coeffs"]
    nb = len(cf)//2
    cfv = -sum(mp.mpf(cf[i])*basis_vals_list[i] + mp.mpf(cf[i+nb])*rt*basis_vals_list[i] for i in range(nb))/a0
    return cfv

# ============ B. Task 1: (sqrt2-1)^2 ============
print("=== B. c0 = (sqrt2-1)^2 ===")
mp.mp.dps = 400
b = mp.sqrt(2)-1; rt2 = mp.sqrt(2); c = b*b
M3, W4, W3, W2 = L7.M3_weights(c)
M3s = mp.mpf(open(WD+"/b2sq_M3_hpl_1150.txt").read())
W4s = mp.mpf(open(WD+"/b2sq_W4_hpl_1150.txt").read())
W3s = mp.mpf(open(WD+"/b2sq_W3_hpl_1150.txt").read())
W2s = mp.mpf(open(WD+"/b2sq_W2_hpl_1150.txt").read())
report("B1 HPL(dps 400) = stored (dps 1150) values",
       abs(M3-M3s) < mp.mpf('1e-380') and abs(W4-W4s) < mp.mpf('1e-380')
       and abs(W3-W3s) < mp.mpf('1e-380') and abs(W2-W2s) < mp.mpf('1e-380'),
       mp.nstr(max(abs(M3-M3s), abs(W4-W4s), abs(W3-W3s), abs(W2-W2s)), 3))
mp.mp.dps = 60
M3q = L7.corner_coeffs_hp(2*mp.acos(c), 3)[3]
report("B2 quadrature(dps 60) vs stored HPL M3", abs(M3q - M3s) < mp.mpf('1e-55'), mp.nstr(abs(M3q - M3s), 3))
mp.mp.dps = 420
W2cf = -Li(2, c) - Li(2, -c) + pi2/12
report("B3 W2 = -Li2(c)-Li2(-c)+pi^2/12 (= -1/2 Li2(b^4)+pi^2/12)", abs(W2s - W2cf) < mp.mpf('1e-400'), mp.nstr(abs(W2s-W2cf), 3))
L2l = mp.log(1+rt2); lam = mp.log(2)
Li3b, Li3mb, Li3s = Li(3, b), Li(3, -b), Li(3, 1/rt2)
Li2b = Li(2, b)
B3v = [z3, Li3b, Li3mb, Li3s, L2l*Li2b, lam*Li2b, pi2*L2l, pi2*lam, lam**3, lam**2*L2l, lam*L2l**2, L2l**3]
W3cf = check_relation("b2sq", "W3", B3v, rt2)
report("B4 W3 stored relation (was identical dps 700/1000)", abs(W3s - W3cf) < mp.mpf('1e-380'), mp.nstr(abs(W3s-W3cf), 3))
pi4 = mp.pi**4
B4v = [Li(4, b), Li(4, -b), Li(4, -b*b), Li(4, 1/rt2), Li(4, mp.mpf(1)/2), Li(4, 1-1/rt2), Li(4, 2-rt2),
       L2l*Li3b, L2l*Li3mb, L2l*Li3s, L2l*z3, lam*Li3b, lam*Li3mb, lam*Li3s, lam*z3,
       Li2b**2, Li2b*L2l**2, Li2b*lam*L2l, Li2b*lam**2, pi2*Li2b,
       pi4, pi2*L2l**2, pi2*lam*L2l, pi2*lam**2,
       L2l**4, L2l**3*lam, L2l**2*lam**2, L2l*lam**3, lam**4]
W4cf = check_relation("b2sq", "W4", B4v, rt2)
report("B5 W4 stored relation (was identical dps 700/1000)", abs(W4s - W4cf) < mp.mpf('1e-380'), mp.nstr(abs(W4s-W4cf), 3))
report("B6 assembly |M3 - (W4cf+W3cf+W2cf)|", abs(M3s - (W4cf+W3cf+W2cf)) < mp.mpf('1e-395'),
       mp.nstr(abs(M3s - (W4cf+W3cf+W2cf)), 3))
# Lie class
QF2 = QuadField(2)
P = QuadPoint(QF2, "(sqrt2-1)^2", (Fr(3), Fr(-2)), primes=[2])
args = {"b": (Fr(-1), Fr(1)), "-b": (Fr(1), Fr(-1)), "-b^2": (Fr(-3), Fr(2)),
        "1/sqrt2": (Fr(0), Fr(1, 2)), "1/2": (Fr(1, 2), Fr(0)),
        "1-1/sqrt2": (Fr(1), Fr(-1, 2)), "2-sqrt2": (Fr(2), Fr(-1))}
coef = {"b": (Fr(-478), Fr(-153)), "-b": (Fr(-470), Fr(-153)), "-b^2": (Fr(83,2), Fr(59,4)),
        "1/sqrt2": (Fr(4), Fr(-4)), "1/2": (Fr(103,2), Fr(41,4)),
        "1-1/sqrt2": (Fr(16), Fr(14)), "2-sqrt2": (Fr(116), Fr(30))}
comp0 = [Fr(0)]*3; comp1 = [Fr(0)]*3
for nm, x in args.items():
    v = li4_coords(P, x)
    p, q = coef[nm]
    for j in range(3):
        comp0[j] += p*v[j]; comp1[j] += q*v[j]
cls = omega_class(P)
report("B7 Li4 part of W4 = class(Omega) = COR7 class (84,-100,128),(-4,-16,48)",
       tuple(comp0) == tuple(cls[0]) and tuple(comp1) == tuple(cls[1])
       and tuple(cls[0]) == (84, -100, 128) and tuple(cls[1]) == (-4, -16, 48))

# ============ C. Task 2: golden corners ============
print("=== C. golden corners ===")
mp.mp.dps = 400
r5 = mp.sqrt(5); phi = (1+r5)/2
L = mp.log(phi); lam = mp.log(2)
for tag, cc in [("gol1", 1/phi), ("gol3", r5-2)]:
    M3, W4, W3, W2 = L7.M3_weights(cc)
    M3s = mp.mpf(open(WD+f"/{tag}_M3_hpl_1150.txt").read())
    W4s = mp.mpf(open(WD+f"/{tag}_W4_hpl_1150.txt").read())
    W3s = mp.mpf(open(WD+f"/{tag}_W3_hpl_1150.txt").read())
    W2s = mp.mpf(open(WD+f"/{tag}_W2_hpl_1150.txt").read())
    report(f"C-{tag}1 HPL(dps 400) = stored",
           abs(M3-M3s) < mp.mpf('1e-380') and abs(W4-W4s) < mp.mpf('1e-380')
           and abs(W3-W3s) < mp.mpf('1e-380') and abs(W2-W2s) < mp.mpf('1e-380'))
    mp.mp.dps = 60
    M3q = L7.corner_coeffs_hp(2*mp.acos(cc), 3)[3]
    report(f"C-{tag}2 quadrature(dps 60) vs stored", abs(M3q - M3s) < mp.mpf('1e-55'), mp.nstr(abs(M3q-M3s), 3))
    mp.mp.dps = 420
    W2cf = -Li(2, cc) - Li(2, -cc) + pi2/12
    report(f"C-{tag}3 W2 = -Li2(c)-Li2(-c)+pi^2/12", abs(W2s - W2cf) < mp.mpf('1e-400'), mp.nstr(abs(W2s-W2cf), 3))
    pi4 = mp.pi**4
    W3cf = check_relation(tag, "W3", [basis_val(nm, phi, r5, L, lam) for nm in json.load(open(WD+f"/{tag}_W3_relation.json"))["basis"]], r5)
    report(f"C-{tag}4 W3 stored relation", abs(W3s - W3cf) < mp.mpf('1e-380'), mp.nstr(abs(W3s-W3cf), 3))
    W4cf = check_relation(tag, "W4", [basis_val(nm, phi, r5, L, lam) for nm in json.load(open(WD+f"/{tag}_W4_relation.json"))["basis"]], r5)
    report(f"C-{tag}5 W4 stored relation", abs(W4s - W4cf) < mp.mpf('1e-400'), mp.nstr(abs(W4s-W4cf), 3))
    report(f"C-{tag}6 assembly", abs(M3s - (W4cf+W3cf+W2cf)) < mp.mpf('1e-395'),
           mp.nstr(abs(M3s - (W4cf+W3cf+W2cf)), 3))
# gol1 W2 elementary
W2s1 = mp.mpf(open(WD+"/gol1_W2_hpl_1150.txt").read())
report("C-gol1-7 W2 = pi^2/20 + L^2/2 (elementary)", abs(W2s1 - (pi2/20 + L**2/2)) < mp.mpf('1e-400'))
# exact weight-3 relations
report("C8 Li3(1/phi^2) = 4 Li3(1/phi) + 4 Li3(-1/phi) (duplication)",
       abs(Li(3, 1/phi**2) - 4*Li(3, 1/phi) - 4*Li(3, -1/phi)) < mp.mpf('1e-400'))
report("C9 6 zeta3 = 30(Li3(1/phi)+Li3(-1/phi)) + pi^2 L - 5 L^3",
       abs(6*z3 - 30*(Li(3, 1/phi)+Li(3, -1/phi)) - pi2*L + 5*L**3) < mp.mpf('1e-400'))
report("C10 Li4(1/phi^6) = 8 Li4(1/phi^3) + 8 Li4(-1/phi^3) (duplication x8)",
       abs(Li(4, 1/phi**6) - 8*Li(4, 1/phi**3) - 8*Li(4, -1/phi**3)) < mp.mpf('1e-400'))
report("C11 Li2 elementary: Li2(1/phi) = pi^2/10 - L^2, Li2(1/phi^2) = pi^2/15 - L^2",
       abs(Li(2, 1/phi) - (pi2/10 - L**2)) < mp.mpf('1e-400') and abs(Li(2, 1/phi**2) - (pi2/15 - L**2)) < mp.mpf('1e-400'))
# Lie classes
QF5 = QuadField(5)
ARGS5 = {"Li4(1/(2phi^2))": (Fr(3,4), Fr(-1,4)), "Li4(1/phi^2)": (Fr(3,2), Fr(-1,2)),
         "Li4(1/2)": (Fr(1,2), Fr(0)), "Li4(2/phi^2)": (Fr(3), Fr(-1)),
         "Li4(phi/2)": (Fr(1,4), Fr(1,4)), "Li4(1/phi^6)": (Fr(9), Fr(-4)), "Li4(-1/phi^3)": (Fr(2), Fr(-1))}
for tag, c0, coefs, exp0, exp1 in [
    ("1/phi", (Fr(-1,2), Fr(1,2)),
     {"Li4(1/(2phi^2))": (Fr(-5,2), Fr(11,6)), "Li4(1/phi^2)": (Fr(113,16), Fr(-187,96)),
      "Li4(1/2)": (Fr(3,10), Fr(19,30)), "Li4(2/phi^2)": (Fr(-7,2), Fr(3,2)),
      "Li4(phi/2)": (Fr(-2), Fr(-2,3)), "Li4(1/phi^6)": (Fr(-29,72), Fr(1,12)),
      "Li4(-1/phi^3)": (Fr(35,18), Fr(7,6))}, (Fr(111,2), 0, -12), (Fr(-51,2), 0, 12)),
    ("1/phi^3", (Fr(-2), Fr(1)),
     {"Li4(1/(2phi^2))": (Fr(1,6), Fr(1,2)), "Li4(1/phi^2)": (Fr(-163,48), Fr(13,16)),
      "Li4(1/2)": (Fr(29,30), Fr(-11,10)), "Li4(2/phi^2)": (Fr(11,2), Fr(1,2)),
      "Li4(phi/2)": (Fr(11,3), Fr(1)), "Li4(1/phi^6)": (Fr(-1,8), Fr(1,72)),
      "Li4(-1/phi^3)": (Fr(1), Fr(7,9))}, (84, -21, 6), (-18, 3, 0))]:
    P = QuadPoint(QF5, tag, c0, primes=[2])
    comp0 = [Fr(0)]*3; comp1 = [Fr(0)]*3
    for nm, (p, q) in coefs.items():
        v = li4_coords(P, ARGS5[nm])
        for j in range(3):
            comp0[j] += p*v[j]; comp1[j] += q*v[j]
    cls = omega_class(P)
    report(f"C12 Lie class {tag} = COR7 printed class",
           tuple(comp0) == tuple(cls[0]) and tuple(comp1) == tuple(cls[1])
           and tuple(cls[0]) == tuple(map(Fr, exp0)) and tuple(cls[1]) == tuple(map(Fr, exp1)))

# ============ D. Task 3: census ============
print("=== D. census ===")
from cor9_rank import true_rank, letter_vecs, rank_int
for D, a, b, e in [(2, 3, -2, 1), (2, -1, 1, 1), (2, 0, 1, 2), (5, -2, 1, 1), (5, -1, 1, 2)]:
    r, *_ = true_rank(D, a, b, e)
    report(f"D1 rank({a},{b})/{e} in Q(sqrt{D}) = 2", r == 2)
for D, a, b, e, rr in [(2, 17, -12, 1, 3), (5, 3, -1, 2, 3), (3, 0, 1, 2, 3), (17, -4, 1, 1, 3), (2, 9, -4, 7, 3)]:
    r, *_ = true_rank(D, a, b, e)
    report(f"D2 rank({a},{b})/{e} in Q(sqrt{D}) = {rr}", r == rr)
# cor7_pt bug example: (2,9,-4,7) correct rank is 3 (buggy letter_span_rank says 4)
vecsM = [[2, 0, 0], [0, -1, 1], [2, -1, 0], [5, -1, 0]]
report("D3 rank of [[2,0,0],[0,-1,1],[2,-1,0],[5,-1,0]] is 3 (cor7_pt buggy elimination says 4)",
       rank_int(vecsM) == 3)
# delta_22 = 0 for the 5 points
QF = {2: QuadField(2), 5: QuadField(5)}
ok = True
for D, a, b, e in [(2, 3, -2, 1), (2, -1, 1, 1), (2, 0, 1, 2), (5, -2, 1, 1), (5, -1, 1, 2)]:
    P = QuadPoint(QF[D], "x", (Fr(a, e), Fr(b, e)), primes=[2])
    d22 = ft_delta22(P.omega_symbol(W4TERMS), 2, len(P.gensE))
    ok = ok and len(d22) == 0
report("D4 true delta_22 = 0 for all 5 rank-2 points", ok)
# census filter sanity: survivors on D=2 match stored
import pickle
surv = pickle.load(open(WD+"/census_filter_survivors_v2.pkl", "rb"))
report("D5 stored survivors (corrected filter): D=2 has 49, D=5 has 28, total 862 over 121 squarefree D<=200",
       len(surv.get(2, [])) == 49 and len(surv.get(5, [])) == 28
       and sum(len(v) for v in surv.values()) == 862)
rank2 = json.load(open(WD+"/census_rank2_complete.json"))
report("D6 census found exactly the 5 rank-2 points",
       sorted(map(tuple, rank2)) == sorted([(2, 3, -2, 1, 2), (2, -1, 1, 1, 2), (2, 0, 1, 2, 2), (5, -2, 1, 1, 2), (5, -1, 1, 2, 2)]))

print(f"\nTOTAL: {sum(1 for _, ok in PASS if ok)}/{len(PASS)} PASS, time {time.time()-T0:.0f}s")
if not all(ok for _, ok in PASS):
    sys.exit(1)
