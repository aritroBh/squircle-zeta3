# 05_reduce.py -- exact reduction of the real dilogarithm values in P1 and P2a to Li2(b), Li2(-b^4)
# by a short explicit chain of named classical identities (real arguments, explicit validity domains):
#   (R) reflection  Li2(x) + Li2(1-x) = pi^2/6 - log x log(1-x)                 0 < x < 1
#   (I) inversion   Li2(x) + Li2(1/x) = -pi^2/6 - (1/2) log^2(-x)               x < 0
#   (L) Landen      Li2(x) + Li2(x/(x-1)) = -(1/2) log^2(1-x)                   x < 1
#   (D) duplication Li2(x) + Li2(-x) = (1/2) Li2(x^2)                           |x| <= 1
#   (A) Abel        Li2(x/(1-y)) + Li2(y/(1-x)) - Li2(xy/((1-x)(1-y))) - Li2(x) - Li2(y)
#                        = log(1-x) log(1-y)                                     x<1, y<1, x+y<1
# Every instance used is printed and checked numerically to 100 digits; the elementary terms are exact
# (logs of S-units of Q(sqrt2) as combinations of lambda=log2, ell3=log3, L=log(1+sqrt2)).
# Then: P1, P2a reduced exactly; P2c from 02; Q0 = 6(P1-P2a-P2c) compared with the theorem exactly.
import sympy as sp, mpmath as mp, json, itertools, sys
from common import lam, l3, L, al, PI, Li2, li2, logS, num, show, Q2, A, B, SQ2, q2, sunit
mp.mp.dps = 110
ref = json.load(open("reference_100.json"))
ld = lambda f: sp.sympify(json.load(open(f))[list(json.load(open(f)).keys())[0]])
P1u, P2au, P2c = ld("P1_unreduced.json"), ld("P2a_unreduced.json"), ld("P2c_exact.json")

def needed_points(*exprs):
    pts = set()
    for e in exprs:
        for f in e.atoms(sp.Function):
            if f.func == Li2: pts.add(q2(f.args[0]))
    return pts
need = needed_points(P1u, P2au)
print("points needed:", sorted(need, key=float))

# ---- pool: x real, x < 1, x and 1-x S-units
IMAX, JMAX, KMAX = [int(t) for t in (sys.argv[1:4] if len(sys.argv) > 3 else (8, 2, 6))]
pool = set()
for sg in (1, -1):
    for i in range(-IMAX, IMAX + 1):
        for j in range(-JMAX, JMAX + 1):
            for k in range(-KMAX, KMAX + 1):
                x = SQ2**i*Q2(3)**j*A**k*sg
                if x < 1 and sunit(1 - x) is not None: pool.add(x)
assert need <= pool, need - pool
print("pool size:", len(pool))
plist = sorted(pool, key=float)

# ---- instances: (name, {point: coeff}, elementary rhs)
inst = []
P2 = PI**2
for x in plist:
    if 0 < x < 1 and (1 - x) in pool:
        inst.append(("R", "x=%r" % x, {x: 1, 1 - x: 1} if x != 1 - x else {x: 2}, P2/6 - logS(x)*logS(1 - x)))
    if x < 0 and (1/x) in pool and x != -1:
        inst.append(("I", "x=%r" % x, {x: 1, 1/x: 1}, -P2/6 - logS(-x)**2/2))
    if x < 1 and x != 0:
        y = x/(x - 1)
        if y in pool and y != x:
            inst.append(("L", "x=%r" % x, {x: 1, y: 1}, -logS(1 - x)**2/2))
    if x != 0 and (x*x) in pool and (-x) in pool and abs(float(x)) <= 1 and x > 0:
        inst.append(("D", "x=%r" % x, {x: 1, -x: 1, x*x: sp.Rational(-1, 2)}, 0))
for x, y in itertools.combinations_with_replacement(plist, 2):
    if x == 0 or y == 0 or not (x + y < 1): continue
    if sunit(1 - x - y) is None: continue
    u1, u2, u3 = x/(1 - y), y/(1 - x), x*y/((1 - x)*(1 - y))
    if not (u1 in pool and u2 in pool and u3 in pool): continue
    cf = {}
    for pt, cc in ((u1, 1), (u2, 1), (u3, -1), (x, -1), (y, -1)):
        cf[pt] = cf.get(pt, 0) + cc
    cf = {p_: c_ for p_, c_ in cf.items() if c_ != 0}
    if cf: inst.append(("A", "x=%r, y=%r" % (x, y), cf, logS(1 - x)*logS(1 - y)))
print("instances:", {t: sum(1 for i_ in inst if i_[0] == t) for t in "RILDA"})

# ---- forward closure: an instance with exactly one unknown point determines it
def closure(basis, goal):
    known = {Q2(0): (0, None), Q2(-1): (-P2/12, "classical Li2(-1) = -pi^2/12"),
             Q2(1, 0)/2: (P2/12 - lam**2/2, "classical Li2(1/2) = pi^2/12 - log^2(2)/2")}
    for bp in basis: known[bp] = (li2(bp), "basis")
    order = {"R": 0, "I": 0, "L": 0, "D": 1, "A": 2}
    changed = True
    while changed and not goal <= set(known):
        changed = False
        for pri in (0, 1, 2):
            for t, lab, cf, rhs in inst:
                if order[t] != pri: continue
                unk = [p_ for p_ in cf if p_ not in known]
                if len(unk) != 1: continue
                u = unk[0]
                val = (rhs - sum(cc*known[p_][0] for p_, cc in cf.items() if p_ != u))/cf[u]
                known[u] = (sp.expand(val), (t, lab, cf, rhs))
                changed = True
            if changed: break
    return known

basis = [B, -(B**4)]
# The two points b^2/sqrt2 and -a^2/sqrt2 (from P2a) are linked to the rest only through one Abel
# instance in which both occur; they enter P2a only through D = Li2(b^2/sqrt2) - Li2(-a^2/sqrt2).
XP, XM = B*B/SQ2, -(A*A)/SQ2
xA, yA = -(A*A), 2*B
abel_star = [i_ for i_ in inst if i_[0] == "A" and i_[1] == "x=%r, y=%r" % (xA, yA)]
assert len(abel_star) == 1
abel_star = abel_star[0]
assert set(abel_star[2]) == {XP, XM, -(A*A), 2*B, A**4*(-1)} and abel_star[2][XP] == 1 and abel_star[2][XM] == -1
need_main = (need - {XP, XM}) | {-(A*A), 2*B, -(A**4)}
known = closure(basis, need_main)
missing = need_main - set(known)
assert not missing, missing
print("basis used: Li2(b), Li2(-b^4)  (no Li2(1/3) needed)")
Dval = sp.expand(abel_star[3] - sum(cc*known[p_][0] for p_, cc in abel_star[2].items() if p_ not in (XP, XM)))
print("Abel* (x=-a^2, y=2b): Li2(b^2/sqrt2) - Li2(-a^2/sqrt2) =", Dval)

# ---- extract the derivation chain actually used (dependencies of the needed points)
used = []
seen = set()
def visit(p_):
    if p_ in seen: return
    seen.add(p_)
    val, why = known[p_]
    if isinstance(why, tuple):
        for q_ in why[2]:
            if q_ != p_: visit(q_)
        used.append((p_, why))
for p_ in sorted(need_main, key=float): visit(p_)
used.append(("D=Li2(b^2/sqrt2)-Li2(-a^2/sqrt2)", abel_star[:4]))
names = {"R": "reflection", "I": "inversion", "L": "Landen", "D": "duplication", "A": "Abel five-term"}
print("\n== derivation chain (%d identity instances); each checked numerically to 100 digits" % len(used))
mx = mp.mpf(0)
for n_, (p_, (t, lab, cf, rhs)) in enumerate(used, 1):
    lhs_num = sum(cc*mp.polylog(2, pt.mp()) for pt, cc in cf.items())
    d_ = abs(lhs_num - num(rhs)); mx = max(mx, d_)
    assert d_ < mp.mpf(10)**-100, (t, lab)
    terms = " + ".join("%s*Li2(%r)" % (cc, pt) for pt, cc in cf.items())
    print("(%2d) %-15s %-40s gives %s" % (n_, names[t], lab, ("Li2(%r)" % p_) if isinstance(p_, Q2) else p_))
    print("      %s = %s      [num. residual %s]" % (terms, rhs, mp.nstr(d_, 2)))
    print("      => %s = %s" % (("Li2(%r)" % p_) if isinstance(p_, Q2) else p_, known[p_][0] if isinstance(p_, Q2) else Dval))
print("max residual of the instances used:", mp.nstr(mx, 3))

# ---- exact reduction of P1, P2a
# P2a: replace Li2(XP) by D + Li2(XM); Li2(XM) must then drop out
P2a_D = sp.expand(P2au.subs(li2(XP), sp.Symbol('Dsym') + li2(XM)))
assert P2a_D.coeff(li2(XM)) == 0
P2a_D = P2a_D.subs(sp.Symbol('Dsym'), Dval)
sub = {li2(p_): known[p_][0] for p_ in need_main if p_ in known}
P1 = sp.expand(P1u.subs(sub)); P2a = sp.expand(P2a_D.subs(sub))
assert not (needed_points(P1, P2a) - {B, -(B**4)}), needed_points(P1, P2a)
r2 = sp.sqrt(2); Lb = li2(B); Lm = li2(-(B**4))
T1 = 3*r2*Lb - lam**2/2 - r2/2*lam*L + (1 + 3*r2/2)*L**2 - (sp.Rational(1, 48) + 3*r2/16)*PI**2
T2a = r2*(3*Lb - Lm/4 - lam*L/2 + L*l3 + L**2/2 - 5*PI**2/24)
T2c = r2*(Lb - Lm/8 + L*l3/2 + PI**2/96 - PI*al/4)
print("\n== exact results")
for nm, val, tg in (("P1", P1, T1), ("P2a", P2a, T2a), ("P2c", P2c, T2c)):
    dd = sp.expand(val - sp.expand(tg))
    print("  %-4s = %s" % (nm, sp.collect(val, [Lb, Lm, PI**2], evaluate=True)))
    print("  %-4s - [closed form of Lemma g3:lem:P]  (exact) = %s" % (nm, dd))
    assert dd == 0
    show("  %s vs 100-digit quadrature" % nm, num(val), mp.mpf(ref[nm]), 95)
Q0 = sp.expand(6*(P1 - P2a - P2c))
Qth = (-PI**2/8 + r2/16*PI**2 + 3*r2/2*PI*al - 3*lam**2 - 9*r2*L*l3 + (6 + 6*r2)*L**2 - 6*r2*Lb + 9*r2/4*Lm)
dq = sp.expand(Q0 - sp.expand(Qth))
print("\n  Q0 = 6(P1-P2a-P2c) =", Q0)
print("  Q0 - [closed form of Theorem thm:S]  (exact) =", dq)
assert dq == 0
mp.mp.dps = 110
Q0n = num(Qth)
print("  Q0 (closed form, 100 digits) =", mp.nstr(Q0n, 100))
show("  Q0 closed form vs 6(P1-P2a-P2c) by quadrature", Q0n, mp.mpf(ref["Q0_from_pieces"]), 95)
for f_, key, tol in (("Q0_k3.json", "Q0_k3", 70), ("Q0_simplex_80.json", "Q0_simplex", 70), ("Q0_simplex_65.json", "Q0_simplex", 60)):
    try:
        vv = mp.mpf(json.load(open(f_))[key])
        show("  Q0 closed form vs %s" % f_, Q0n, vv, tol)
    except FileNotFoundError:
        print("  (%s not present)" % f_)
json.dump({"chain": [(repr(p_), t, lab) for p_, (t, lab, cf, rhs) in used],
           "P1": str(P1), "P2a": str(P2a), "P2c": str(P2c), "Q0": str(Q0)}, open("05_reduce_out.json", "w"), indent=1)
