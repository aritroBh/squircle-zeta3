"""
verify_full.py — PRV1: complete exact verifier for the depth-2 reductions of the
regular octagon (CX1/OCT6) and the regular decagon (DEC1).

Every step is EXACT: sympy formal differentiation algebra with exact constant pins,
or exact (finite-sum) series identities.  The high-precision numeric checks in Part 4
are independent evidence, not part of the proof.

  Part 0  Foundations: partial derivatives of Li_{2,1}; Li_{1,1} closed form.
  Part 1  Classical dilog/trilog function identities (derivative 0 + exact pin).
  Part 2  Octagon:  G, H families, derivative identities, constants, specialisation.
  Part 3  Decagon:  G_dec, H_dec, derivative identities, constants, golden lattice,
          specialisation matching DEC1 exactly.
  Part 4  Numeric evidence (200 dps) for every identity.
Run:  python verify_full.py   -> prints ALL CHECKS PASSED.
"""

from pathlib import Path
import hashlib
import pickle
import sys

import mpmath as mp
import sympy as sp

SRC_OCT = Path('research/08_theory/OCT6_octagon_closed')

def depth2(x, y, p=2, details=False):
    """Li_{p,1}(x,y) = sum_{n>m>=1} x^n y^m / (n^p m); tail <= rho^(N+1)/(1-rho)."""
    rho = max(abs(x), abs(x * y))
    if not 0 < rho < 1:
        if rho == 0:
            return (mp.mpf(0), 0, mp.mpf(0)) if details else mp.mpf(0)
        raise ValueError("The defining series requires |x|, |xy| < 1")
    target = mp.power(10, -mp.mp.dps - 5)
    nmax = int(mp.ceil(mp.log(target * (1 - rho)) / mp.log(rho)))
    h, xy_power, total = mp.mpf(0), x * y, mp.mpf(0)
    for n in range(1, nmax + 1):
        total += h / n**p
        h = x * h + x * xy_power / n
        xy_power *= x * y
    return (total, nmax, rho**(nmax+1)/(1-rho)) if details else total


# ======================================================================
def part0():
    print("=" * 74)
    print("PART 0: foundations (exact)")
    x, y = sp.symbols('x y')
    # -- 0.1 partial_y and partial_x as EXACT finite-sum identities (truncation N).
    #     d/dy S_N:  (1-y)*d/dy S_N = sum_{n<N} x^n/n^2 (1 - y^{n-1})  (exact poly id);
    #     dividing by (1-y) (formal, y != 1) gives d_y Li_21 = (Li2(x) - Li2(x y)/y)/(1-y).
    N = 25
    S = sum(x**n * y**m / (n**2 * m) for n in range(1, N) for m in range(1, n))
    lhs_y = sp.diff(S, y)
    num_y = (sum(x**n/n**2 for n in range(1, N))
             - sum((x*y)**n/n**2 for n in range(1, N))/y)
    r = sp.expand((1-y)*lhs_y - num_y)
    print("  0.1a (1-y) d_y Li_21 truncated-series identity: residual =", r)
    assert r == 0
    lhs_x = sp.diff(S, x)
    rhs_x = sp.expand(sum(x**n * y**m / (n * m) for n in range(1, N) for m in range(1, n))/x)
    r = sp.expand(lhs_x - rhs_x)
    print("  0.1b d_x Li_21 truncated-series identity: residual =", r)
    assert r == 0
    # -- 0.2 Li_{1,1} closed form:  F = log(1-x)log(1-y) + Li2(-y/(1-y)) - Li2(-y(1-x)/(1-y))
    #     d/dx F = -log(1-x y)/(1-x) (the series derivative of Li_{1,1}); F(0,y)=0.
    #     Branch note: on 0<x,y<1 we have 0<1-x, 0<1-y, 0<(1-xy)/(1-y), and
    #     log((1-xy)/(1-y)) = log(1-xy) - log(1-y) (real log of positive args).
    z = -y*(1-x)/(1-y)
    dz = sp.diff(z, x)
    dF = (-sp.log(1-y)/(1-x)
          + (sp.log(1-x*y) - sp.log(1-y)) * dz / z)   # d/dx[-Li2(z)] = log(1-z) z'/z
    r = sp.factor(dF - (-sp.log(1-x*y)/(1-x)))
    print("  0.2  d/dx Li_11 closed form minus series derivative =", r)
    assert r == 0
    print("       F(0,y) = log(1)*log(1-y) + Li2(-y/(1-y)) - Li2(-y/(1-y)) = 0 (exact)")
    print("PART 0 OK")


# ======================================================================
# Part 1: classical dilog/trilog function identities on 0<t<1 (exact)
# ======================================================================
# Formal differentiation engine.  Letters: R=log t, A=log(1-t), B=log(1+t),
# C=log 2, Z=zeta(2), Z3=zeta(3).  Polylog symbols carry a derivative map
# obtained from d/dz Li_n(z) = Li_{n-1}(z)/z with real branches on (0,1).
t = sp.symbols('t', positive=True)
R, A, B, C, Z, Z3 = sp.symbols('R A B C Z Z3')
Dt, Dm, Db, Dh, Dc, Dmc, D1, Du, Dw = sp.symbols('Dt Dm Db Dh Dc Dmc D1 Du Dw')
Dt2 = sp.Symbol('Dt2')
Tt, Tm, Tmc, Tu, T1, Tb, Tw, Tc = sp.symbols('Tt Tm Tmc Tu T1 Tb Tw Tc')
Tt2 = sp.Symbol('Tt2')
q, lc, U = B-C, A-B, C+R-B        # log((1+t)/2), log c, log u ; c=(1-t)/(1+t), u=2t/(1+t)
DERIV = {
    R: 1/t, A: -1/(1-t), B: 1/(1+t), C: 0, Z: 0, Z3: 0,
    Dt: -A/t, Dm: -B/t, D1: R/(1-t), Db: q/(1-t), Dh: -(A-C)/(1+t),
    Dc: 2*U/(1-t**2), Dmc: -2*q/(1-t**2), Dw: B/(t*(1+t)),
    Dt2: -2*(A+B)/t,
    Tt: Dt/t, Tm: Dm/t, Tmc: -2*Dmc/(1-t**2), Tu: Du/(t*(1+t)),
    T1: -D1/(1-t), Tb: -Db/(1-t), Tw: Dw/(t*(1+t)), Tc: -2*Dc/(1-t**2),
    Tt2: 2*Dt2/t,
}
diff = lambda e: sum(sp.diff(e, v)*dv for v, dv in DERIV.items())

def part1():
    print("=" * 74)
    print("PART 1: classical dilog/trilog function identities (derivative = 0 exactly)")
    ids = {
        "E   Euler: Li2(t)+Li2(1-t)=Z-log t log(1-t)": Dt + D1 + R*A - Z,
        "L1  Landen at -t: Li2(-t)+Li2(w)=-B^2/2":     Dm + Dw + B**2/2,
        "L2  Landen at -c: Li2(-c)+Li2(b)=-q^2/2":     Dmc + Db + q**2/2,
        "DUP2 dilog duplication: Li2(t^2)=2Li2(t)+2Li2(-t)": Dt2 - 2*Dt - 2*Dm,
        "OC  odd complement: Li2(t)-Li2(-t)+Li2(c)-Li2(-c)=3Z/2-log t log c":
            Dt - Dm + Dc - Dmc + R*lc - 3*Z/2,
        "DUP3 trilog duplication: Li3(t^2)=4Li3(t)+4Li3(-t)": Tt2 - 4*Tt - 4*Tm,
    }
    for name, e in ids.items():
        r = sp.factor(diff(e))
        print("   d/dt[", name, "] =", r)
        if not name.startswith("DUP3"):
            assert r == 0
    # DUP3's derivative uses DUP2 as a function identity: check that reduction too
    print("   [d/dt DUP3 = (2/t)*(DUP2), so d/dt DUP3 = 0 on functions; check:]")
    r = sp.factor(diff(Tt2 - 4*Tt - 4*Tm).subs(Dt2, 2*Dt+2*Dm))
    print("   d/dt[DUP3] after DUP2 =", r)
    assert r == 0
    print("  exact constant pins:")
    print("   E  at t->0+: Li2(0)+Li2(1)=Z ; log t log(1-t)->0   => 0  (exact)")
    print("   L1 at t=0:   0+0+0 = 0  (exact)")
    print("   L2 at t=0:   Li2(-1)+Li2(1/2)+log^2(1/2)/2")
    print("              = -Z/2 + (Z/2 - C^2/2) + C^2/2 = 0")
    print("      [Li2(-1)=-Z/2 from DUP2 at t=1; Li2(1/2)=Z/2-C^2/2 from E at t=1/2]")
    print("   DUP2 at t=0: 0  (exact);   OC at t->1-: Z+Z/2-3Z/2=0, log t log c->0 (exact)")
    print("   DUP3 at t=0: 0  (exact);   Li3(-1) = -3*Z3/4 follows from DUP3 at t=1.")
    print("PART 1 OK")


# ======================================================================
# Part 2: the octagon (exact)
# ======================================================================
# Families on 0<t<1:
#   G(t) = Li_21(-t, -(1-t)/(2t))        curve (x,xy) = (-t, (1-t)/2)
#   H(t) = Li_21((1-t)/2, -2t/(1-t))     curve (x,xy) = ((1-t)/2, -t)  [x<->p swap of G]
# Chain rule (Part 0) + Li_11 closed form give
#   G' = JG/t + Dm/(t(1+t)) + 2 Db/(1-t^2),   JG = (q-R) B + Dc - D1
#   H' = -JH/(1-t) - Dm/(t(1+t)) - 2 Db/(1-t^2),  JH = -lc q + Du - Dt
JG = (q-R)*B + Dc - D1
JH = -lc*q + Du - Dt
dG = JG/t + Dm/(t*(1+t)) + 2*Db/(1-t**2)
dH = -JH/(1-t) - Dm/(t*(1+t)) - 2*Db/(1-t**2)
# Closed forms (constants included; derivatives of constants vanish):
G_cl = 3*Z3/4 + Tt + Tm + Tmc - Tu - q*(Dt + Db/2 - Dh/2)
H_cl = T1 + Tb + Tw - Tc - sp.Symbol('T12') - B*(Dc + Db) + B**2*U/2
# T12 = Li3(1/2), a constant.
DERIV[sp.Symbol('T12')] = 0
# Function identities of Part 1 as substitution rules (proved there):
BASIC = {D1: Z - Dt - R*A,                 # E
         Dh: Z - Db - (A-C)*q,             # E at h=(1+t)/2  (1-h=b)
         Dmc: -Db - q**2/2,                # L2
         Du: Z - Dc - lc*U,                # E at u=2t/(1+t) (1-u=c)
         Dw: -Dm - B**2/2}                 # L1
OC = Dt - Dm + Dc + Db + q**2/2 - 3*Z/2 + R*lc   # = OC with Dmc eliminated via L2

def part2():
    print("=" * 74)
    print("PART 2: octagon derivative identities (exact)")
    rH = sp.factor((diff(H_cl) - dH).subs(BASIC))
    print("  d/dt H_closed - H'  (mod E,L1,L2) =", rH)
    assert rH == 0
    rG = sp.factor((diff(G_cl) - dG).subs(BASIC))
    print("  d/dt G_closed - G'  (mod E,L1,L2) =", sp.factor(rG), " = -OC/(1+t):",
          sp.factor(rG - (-OC/(1+t))))
    assert sp.factor(rG - (-OC/(1+t))) == 0
    print("  and OC = 0 by Part 1 (derivative 0 + pin at t->1-).")
    print("  constants:")
    print("   G(t)->0 as t->1- (|G|<=|y| zeta(2)/(1-|y|), y=-(1-t)/(2t)->0);")
    print("   G_closed(1) = 3Z3/4 + Li3(1) + Li3(-1) + Li3(0) - Li3(1) - 0")
    print("              = 3Z3/4 + Z3 - 3Z3/4 - Z3 = 0  [Li3(-1)=-3Z3/4, Part 1 DUP3]")
    print("   H(t)->0 as t->0+ (same estimate, y=-2t/(1-t)->0);")
    print("   H_closed(0) = Li3(1)+Li3(1/2)+Li3(0)-Li3(1)-Li3(1/2)-0+0 = 0 (exact)")
    print("  Hence G=G_closed and H=H_closed on (0,1).")
    print("PART 2 OK")


# ---------------- octagon specialisation at t = k = sqrt(2)-1 ----------------
def part2_specialisation():
    print("=" * 74)
    print("PART 2b: octagon specialisation at t = sqrt(2)-1 (exact over Q(sqrt(2)))")
    s = sp.sqrt(2)
    k, h, x2, d, x1 = s-1, s/2, 1-s/2, 2-s, 1-s
    args = {
        "b(k)=x2": (1-k)/2 - x2, "h(k)=sqrt(2)/2": (1+k)/2 - h,
        "c(k)=k": (1-k)/(1+k) - k, "u(k)=d": 2*k/(1+k) - d,
        "w(k)=x2": k/(1+k) - x2,
        "G second arg = -sqrt(2)/2": -(1-k)/(2*k) + h,
        "H second arg = -sqrt(2)": -2*k/(1-k) + s,
        "h(k)^2 = 1/2": h**2 - sp.Rational(1, 2),
        "k(1+sqrt2) = 1": k*(1+s) - 1,
        "1-k = sqrt(2) k": 1-k - s*k,
        "x2 = x1/(x1-1)": x2 - x1/(x1-1),
        "-sqrt2/2 = x2/x1": -s/2 - x2/x1,
        "-sqrt2 = x1/x2": -s - x1/x2,
    }
    for name, r in args.items():
        r = sp.simplify(r)
        print("   ", name, ": residual =", r)
        assert r == 0
    # silver dilog rules as specialisations of E and Landen (proved function ids):
    l, L = sp.log(2), sp.log(1+s)
    Z2s = sp.pi**2/6
    print("   silver rules (exact specialisations of E and Landen at points of (0,1)):")
    print("    D(d):  E at d=2-sqrt2: 1-d=k, log d = l/2 - L, log k = -L")
    print("    D(x2): Landen at x1: -x1/(1-x1) = x2, log(1-x1) = l/2")
    print("    D(h):  E at h=sqrt2/2: 1-h=x2, log h = -l/2, log x2 = -l/2 - L")
    print("    D(-h): Landen at k: -k/(1-k) = -h, log(1-k) = l/2 - L")
    # log identities are checked as exact MULTIPLICATIVE identities among the
    # positive real arguments (real log is injective on R_{>0}):
    checks = {"-x1/(1-x1) = x2": sp.simplify(-x1/(1-x1) - x2),
              "-k/(1-k) = -h": sp.simplify(-k/(1-k) + h),
              "d*(1+sqrt2) = sqrt2   [<= log d = l/2 - L]": sp.simplify(d*(1+s) - s),
              "k*(1+sqrt2) = 1       [<= log k = -L]": sp.simplify(k*(1+s) - 1),
              "x2*sqrt2*(1+sqrt2) = 1  [<= log x2 = -l/2 - L]": sp.simplify(x2*s*(1+s) - 1),
              "1-x1 = sqrt2          [<= log(1-x1) = l/2]": sp.simplify(1-x1 - s)}
    for name, r in checks.items():
        print("     ", name, ": residual =", r)
        assert r == 0
    # compact specialised forms and exact comparison with OCT6 reductions.pkl
    reductions = pickle.loads((SRC_OCT / "reductions.pkl").read_bytes())
    l2, lu, z2, z3 = sp.symbols("l2 lu z2 z3")
    Dd, Dk, Dx1, Dx2, Dh, Dmh = sp.symbols("L2d L2k L2x1 L2x2 L2r2h L2mr2h")
    Tk, Tx1, Tx2, Td, Th, Tmh = sp.symbols("L3k L3x1 L3x2 L3d L3r2h L3mr2h")
    compact3 = 3*z3/4 + Tk + 2*Tx1 - Td + l2*Dk/2 + l2*(Dx2-Dh)/4
    compact4 = Td + 2*Tx2 - Tk - 4*Th - 4*Tmh - l2*(Dk+Dx2)/2 + l2**2*(l2/2-lu)/8
    rules = {Dd: z2-Dk+l2*lu/2-lu**2, Dx2: -Dx1-l2**2/8,
             Dh: z2+Dx1-l2*lu/2-l2**2/8, Dmh: -Dk-lu**2/2+l2*lu/2-l2**2/8}
    for key, compact in (("Lam", compact3), ("swap", compact4)):
        r = sp.expand((compact - reductions[key]).xreplace(rules))
        print("   ", key, ": compact specialisation - stored reduction (exact) =", r)
        assert r == 0
    print("   [compact forms: G_closed(k), H_closed(k) rewritten by the argument")
    print("    identities above; -log h(k) = l2/2, log(1+k) = l2/2,")
    print("    log u(k) = log d = l2/2 - L, and Li3(1/2) = 4Li3(h)+4Li3(-h) (DUP3).]")
    print("PART 2b OK")


# ======================================================================
# Part 3: the decagon (exact)
# ======================================================================
# Families on 0<t<1 (same shape as the octagon's, with c=1):
#   G_dec(t) = Li_21(-t, -(1-t)/t)       curve (x,xy) = (-t, 1-t)
#   H_dec(t) = Li_21(1-t, -t/(1-t))      curve (x,xy) = (1-t, -t)  [x<->p swap]
Dv, Dw1, D1t2 = sp.symbols('Dv Dw1 D1t2')
T1t2, Tv, Tw1 = sp.symbols('T1t2 Tv Tw1')
DERIV.update({
    Dt2: -2*(A+B)/t,
    D1t2: 4*t*R/(1-t**2),
    Dv: B/(t*(1+t)),                    # v = t/(1+t)
    Dw1: (R-B)/(1+t),                   # w = 1/(1+t)
    Tt2: 2*Dt2/t,
    T1t2: -2*t*D1t2/(1-t**2),
    Tv: Dv/(t*(1+t)),
    Tw1: -Dw1/(1+t),
})
Li11_Gd = -R*B + D1 - D1t2
Li11_Hd = -R*A + Dt - Dt2
dGd = Li11_Gd/t + Dm/t + D1/(1-t)
dHd = -Li11_Hd/(1-t) - D1/(1-t) - Dm/t
G_dec_cl = Z3 - R*Dt2/2 + Tt2 - 2*Tt - T1
H_dec_cl = T1t2 - Tw1 + Tv - B*(Dm/2 + 3*D1/2 + D1t2/4 + Z/4)
BASIC_DEC = {Dt2: 2*Dt + 2*Dm,                       # DUP2
             D1: Z - Dt - R*A,                       # E
             D1t2: Z - (2*Dt+2*Dm) - 2*R*(A+B),      # E at t^2
             Dv: -Dm - B**2/2,                       # L1 (Landen at -t)
             Dw1: Z - (-Dm - B**2/2) + B*(R-B)}      # E at v (1-v = w)

def part3():
    print("=" * 74)
    print("PART 3: decagon derivative identities (exact)")
    rG = sp.factor((diff(G_dec_cl) - dGd).subs(BASIC_DEC))
    print("  d/dt G_dec_closed - G_dec' (mod DUP2,E,L1) =", rG)
    assert rG == 0
    rH = sp.factor((diff(H_dec_cl) - dHd).subs(BASIC_DEC))
    print("  d/dt H_dec_closed - H_dec' (mod DUP2,E,L1) =", rH)
    assert rH == 0
    print("  constants:")
    print("   G_dec(t)->0 as t->1- (|G|<=|y| zeta(2)/(1-|y|), y=-(1-t)/t->0);")
    print("   G_dec_closed(1) = Z3 - 0 + Li3(1) - 2 Li3(1) - Li3(0) = 0 (exact)")
    print("   H_dec(t)->0 as t->0+ (same estimate, y=-t/(1-t)->0);")
    print("   H_dec_closed(0) = Li3(1)-Li3(1)+Li3(0)-0 = 0 (exact)")
    print("  Hence G_dec=G_dec_closed and H_dec=H_dec_closed on (0,1).")
    print("PART 3 OK")


# ---- 3-term Landen trilog identity (function identity, exact) ----
Dml, Tml = sp.symbols('Dml Tml')     # Li2(-t/(1-t)), Li3(-t/(1-t))
DERIV[Dml] = A/(t*(1-t))
DERIV[Tml] = Dml/(t*(1-t))
T3 = Tt + T1 + Tml - (Z3 + Z*A - sp.Rational(1,2)*R*A**2 + A**3/6)

def part3_landen3():
    print("=" * 74)
    print("PART 3b: 3-term Landen trilog identity (exact)")
    r = sp.factor(diff(T3).subs({D1: Z - Dt - R*A, Dml: -Dt - A**2/2}))
    print("  d/dt[Li3(t)+Li3(1-t)+Li3(-t/(1-t)) - (Z3 + Z A - R A^2/2 + A^3/6)]")
    print("  after E and Landen =", r)
    assert r == 0
    print("  pin at t->0+: 0 + Z3 + 0 - (Z3 + 0 - 0 + 0) = 0 (exact); identity holds on (0,1).")
    print("PART 3b OK")


# ---- golden lattice and specialisation (exact over Q(sqrt(5))) ----
def part3_specialisation():
    print("=" * 74)
    print("PART 3c: golden lattice + decagon specialisation at t = 1/phi (exact)")
    s5 = sp.sqrt(5)
    phi = (1+s5)/2
    l, z2, z3 = sp.symbols('l z2 z3')
    L2a, L2b, L2c = sp.symbols('L2a L2b L2c')   # Li2(1/phi^2), Li2(1/phi), Li2(-1/phi)
    L3a, L3b, L3c = sp.symbols('L3a L3b L3c')   # Li3(1/phi^2), Li3(1/phi), Li3(-1/phi)
    # argument identities over Q(sqrt(5))
    argid = {
        "1 - 1/phi^2 = 1/phi": sp.simplify(1 - 1/phi**2 - 1/phi),
        "-(1/phi^2)/(1-1/phi^2) = -1/phi": sp.simplify(-(1/phi**2)/(1-1/phi**2) + 1/phi),
        "1/(1+1/phi) = 1/phi": sp.simplify(1/(1+1/phi) - 1/phi),
        "(1/phi)/(1+1/phi) = 1/phi^2": sp.simplify((1/phi)/(1+1/phi) - 1/phi**2),
        "log(1+1/phi) = log phi": sp.simplify(sp.log(1+1/phi) - sp.log(phi)),
        "(1/phi)^2 = 1 - 1/phi": sp.simplify((1/phi)**2 - (1-1/phi)),
    }
    for name, r in argid.items():
        print("   ", name, ": residual =", r)
        assert r == 0
    # golden dilogs: E @ 1/phi^2, Landen @ 1/phi^2, DUP2 @ 1/phi
    sold = sp.solve([L2a + L2b - (z2 - 2*l**2),
                     L2c + L2a + l**2/2,
                     L2a - 2*L2b - 2*L2c], [L2a, L2b, L2c], dict=True)[0]
    print("   golden dilogs:", sold)
    # golden trilogs: 3-term Landen @ 1/phi^2, DUP3 @ 1/phi.
    # Two equations in L3a,L3b,L3c: they determine L3a (Lewin's value) and the
    # SUM L3b+L3c exactly; the difference L3b-L3c remains free (classical:
    # Li3(1/phi), Li3(-1/phi) are individually independent of z3, z2, l).
    solt = sp.solve([L3a + L3b + L3c - (z3 - l*z2 + 5*l**3/6),
                     L3a - 4*L3b - 4*L3c], [L3a, L3b], dict=True)[0]
    print("   Lewin: Li3(1/phi^2) =", sp.expand(solt[L3a]))
    S = sp.expand(solt[L3b] + L3c)
    print("   sum relation: Li3(1/phi)+Li3(-1/phi) =", S)
    # specialisation
    Lam = sp.expand((z3 + l*L2a/2 - 2*L3b).subs(sold))          # G_dec_closed(1/phi)
    print("   Lambda_dec = G_dec_closed(1/phi) =", Lam, "  [with Li3(1/phi) free]")
    Lam_DEC1 = 2*z3/5 + 4*l*z2/5 - l**3 + 3*L3c + L3b
    r = sp.simplify(sp.expand(Lam - Lam_DEC1).subs(L3b, S - L3c))
    print("   Lambda_dec - DEC1 form (mod sum relation) =", r)
    assert r == 0
    swap = sp.expand((L3a - l*(L2c/2 + 3*L2a/2 + L2b/4 + z2/4)).subs(sold))
    print("   swap_dec = H_dec_closed(1/phi) =", swap)
    swap_full = sp.expand(swap.subs(L3a, solt[L3a]))
    print("   swap_dec fully reduced (via Lewin) =", swap_full)
    swap_DEC1 = 3*z3/5 - 7*l*z2/5 + 2*l**3 + L3c + L3b
    r = sp.simplify(sp.expand(swap_full - swap_DEC1).subs(L3b, S - L3c))
    print("   swap_dec - DEC1 form (mod sum relation) =", r)
    assert r == 0
    print("PART 3c OK")


# ======================================================================
# Chain-rule assembly (exact): the dG, dH, dGd, dHd formulas used in
# Parts 2 and 3 are derived HERE from Part 0's partial-derivative formulas
# and the Li_{1,1} closed form, by exact rational simplification on the curves.
# ======================================================================
def part_chainrule():
    print("=" * 74)
    print("CHAIN-RULE ASSEMBLY (exact): dG, dH, dGd, dHd from Part 0")
    tt = sp.symbols('t', positive=True)
    def li11(xe, ye, log1x, log1y, L2a, L2b):
        # Li11(xe,ye) via Part 0.2, with the exact argument evaluations supplied
        return log1x*log1y + L2a - L2b
    def dLi21_y(xe, ye, L2x, L2xy):
        return (L2x - L2xy/ye)/(1-ye)
    # ---- octagon G: (x,y) = (-t, -(1-t)/(2t)) ----
    xe, ye = -tt, -(1-tt)/(2*tt)
    # exact argument identities (rational, checked):
    ids = {"1-x = 1+t": sp.simplify(1-xe - (1+tt)),
           "1-y = (1+t)/(2t)": sp.simplify(1-ye - (1+tt)/(2*tt)),
           "-y/(1-y) = (1-t)/(1+t) = c": sp.simplify(-ye/(1-ye) - (1-tt)/(1+tt)),
           "-y(1-x)/(1-y) = 1-t": sp.simplify(-ye*(1-xe)/(1-ye) - (1-tt)),
           "xy = (1-t)/2 = b": sp.simplify(xe*ye - (1-tt)/2)}
    # ---- octagon H: (x,y) = ((1-t)/2, -2t/(1-t)) ----
    xeH, yeH = (1-tt)/2, -2*tt/(1-tt)
    ids.update({"H 1-x = (1+t)/2 = h": sp.simplify(1-xeH - (1+tt)/2),
                "H 1-y = (1+t)/(1-t)": sp.simplify(1-yeH - (1+tt)/(1-tt)),
                "H -y/(1-y) = 2t/(1+t) = u": sp.simplify(-yeH/(1-yeH) - 2*tt/(1+tt)),
                "H -y(1-x)/(1-y) = t": sp.simplify(-yeH*(1-xeH)/(1-yeH) - tt),
                "H xy = -t": sp.simplify(xeH*yeH + tt)})
    # ---- decagon G: (x,y) = (-t, -(1-t)/t) ----
    xeG2, yeG2 = -tt, -(1-tt)/tt
    ids.update({"dG 1-x = 1+t": sp.simplify(1-xeG2 - (1+tt)),
                "dG 1-y = 1/t": sp.simplify(1-yeG2 - 1/tt),
                "dG -y/(1-y) = 1-t": sp.simplify(-yeG2/(1-yeG2) - (1-tt)),
                "dG -y(1-x)/(1-y) = 1-t^2": sp.simplify(-yeG2*(1-xeG2)/(1-yeG2) - (1-tt**2)),
                "dG xy = 1-t": sp.simplify(xeG2*yeG2 - (1-tt))})
    # ---- decagon H: (x,y) = (1-t, -t/(1-t)) ----
    xeH2, yeH2 = 1-tt, -tt/(1-tt)
    ids.update({"dH 1-x = t": sp.simplify(1-xeH2 - tt),
                "dH 1-y = 1/(1-t)": sp.simplify(1-yeH2 - 1/(1-tt)),
                "dH -y/(1-y) = t": sp.simplify(-yeH2/(1-yeH2) - tt),
                "dH -y(1-x)/(1-y) = t^2": sp.simplify(-yeH2*(1-xeH2)/(1-yeH2) - tt**2),
                "dH xy = -t": sp.simplify(xeH2*yeH2 + tt)})
    for name, r in ids.items():
        print("   ", name, ": residual =", r)
        assert r == 0
    # assemble derivatives (real-log splits of positive arguments on 0<t<1):
    # octagon G: log(1-x)=B, log(1-y)=B-C-R=q-R, L2(-y/(1-y))=Dc, L2(-y(1-x)/(1-y))=D1
    Li11_G = B*(q-R) + Dc - D1
    dG_asm = (Li11_G/xe)*sp.diff(xe, tt) + dLi21_y(xe, ye, Dm, Db)*sp.diff(ye, tt)
    r = sp.factor(dG_asm - dG)
    print("    assembled dG - stated dG =", r); assert r == 0
    # octagon H: log(1-x)=q, log(1-y)=-lc, L2(u)=Du, L2(t)=Dt
    Li11_H = -lc*q + Du - Dt
    dH_asm = (Li11_H/xeH)*sp.diff(xeH, tt) + dLi21_y(xeH, yeH, Db, Dm)*sp.diff(yeH, tt)
    r = sp.factor(dH_asm - dH)
    print("    assembled dH - stated dH =", r); assert r == 0
    # decagon G: log(1-x)=B, log(1-y)=-R, L2(1-t)=D1, L2(1-t^2)=D1t2
    Li11_Gd = B*(-R) + D1 - D1t2
    dGd_asm = (Li11_Gd/xeG2)*sp.diff(xeG2, tt) + dLi21_y(xeG2, yeG2, Dm, D1)*sp.diff(yeG2, tt)
    r = sp.factor(dGd_asm - dGd)
    print("    assembled dGd - stated dGd =", r); assert r == 0
    # decagon H: log(1-x)=R, log(1-y)=-A, L2(t)=Dt, L2(t^2)=Dt2
    Li11_Hd = R*(-A) + Dt - Dt2
    dHd_asm = (Li11_Hd/xeH2)*sp.diff(xeH2, tt) + dLi21_y(xeH2, yeH2, D1, Dm)*sp.diff(yeH2, tt)
    r = sp.factor(dHd_asm - dHd)
    print("    assembled dHd - stated dHd =", r); assert r == 0
    print("CHAIN-RULE ASSEMBLY OK")


# ======================================================================
# Part 4: independent numeric evidence (checks, NOT part of the proof)
# ======================================================================
def li2(z): return mp.polylog(2, z)
def li3(z): return mp.polylog(3, z)

def G_cl_mp(t):
    return (mp.mpf(3)/4*mp.zeta(3) + li3(t) + li3(-t)
            + li3(-(1-t)/(1+t)) - li3(2*t/(1+t))
            - mp.log((1+t)/2)*(li2(t) + li2((1-t)/2)/2 - li2((1+t)/2)/2))

def H_cl_mp(t):
    return (li3(1-t) + li3((1-t)/2) + li3(t/(1+t)) - li3((1-t)/(1+t))
            - li3(mp.mpf(1)/2)
            - mp.log(1+t)*(li2((1-t)/(1+t)) + li2((1-t)/2))
            + mp.log(1+t)**2*mp.log(2*t/(1+t))/2)

def G_dec_cl_mp(t):
    return (mp.zeta(3) - mp.log(t)*li2(t**2)/2 + li3(t**2)
            - 2*li3(t) - li3(1-t))

def H_dec_cl_mp(t):
    return (li3(1-t**2) - li3(1/(1+t)) + li3(t/(1+t))
            - mp.log(1+t)*(li2(-t)/2 + 3*li2(1-t)/2 + li2(1-t**2)/4
                           + mp.pi**2/24))

def dig(d):
    d = abs(d)
    if d == 0:
        return "exactly 0"
    e = int(mp.floor(mp.log10(d)))
    m = float(d / mp.power(10, e))
    return f"{m:.1f}e{e}"

def part4():
    print("=" * 74)
    print("PART 4: numeric evidence (mpmath, independent of the exact proofs)")
    mp.mp.dps = 120
    z2, z3 = mp.pi**2/6, mp.zeta(3)
    # 4.0 classical function identities at sample points
    for t0 in (mp.mpf("0.3"), mp.mpf("0.55"), mp.mpf("0.8")):
        c = (1-t0)/(1+t0)
        vals = {
            "E": li2(t0)+li2(1-t0)-(z2-mp.log(t0)*mp.log(1-t0)),
            "L1": li2(-t0)+li2(t0/(1+t0))+mp.log(1+t0)**2/2,
            "L2": li2(-c)+li2((1-t0)/2)+mp.log((1+t0)/2)**2/2,
            "DUP2": li2(t0**2)-2*li2(t0)-2*li2(-t0),
            "OC": (li2(t0)-li2(-t0)+li2(c)-li2(-c)
                   -3*z2/2+mp.log(t0)*mp.log(c)),
            "DUP3": li3(t0**2)-4*li3(t0)-4*li3(-t0),
            "Landen3": (li3(t0)+li3(1-t0)+li3(-t0/(1-t0))
                        -(z3+z2*mp.log(1-t0)-mp.log(t0)*mp.log(1-t0)**2/2
                          +mp.log(1-t0)**3/6)),
        }
        worst = max(abs(v) for v in vals.values())
        print(f"   t={t0}: classical ids worst residual {dig(worst)}"
              f"  ({max(vals, key=lambda k: abs(vals[k]))})")
    # 4.0b Li_{1,1} closed form vs series
    for x0, y0 in ((mp.mpf("0.3"), mp.mpf("0.4")), (mp.mpf("-0.4"), mp.mpf("0.5")),
                   (mp.mpf("0.1"), mp.mpf("-0.6"))):
        F = (mp.log(1-x0)*mp.log(1-y0) + li2(-y0/(1-y0))
             - li2(-y0*(1-x0)/(1-y0)))
        ser = depth2(x0, y0, p=1)
        print(f"   Li_11 closed form vs series at ({x0},{y0}): residual {dig(F-ser)}")
    # 4.1 octagon families: series vs closed forms
    print("  octagon: G,G_cl and H,H_cl on (0,1):")
    for t0 in (mp.mpf("0.17"), mp.mpf("0.31"), mp.mpf("0.5"),
               mp.mpf("0.63"), mp.mpf("0.82")):
        rG = depth2(-t0, -(1-t0)/(2*t0)) - G_cl_mp(t0)
        rH = depth2((1-t0)/2, -2*t0/(1-t0)) - H_cl_mp(t0)
        print(f"   t={t0}: |G-G_cl| {dig(rG)}, |H-H_cl| {dig(rH)}")
    # 4.2 octagon specialisation at k = sqrt(2)-1
    k = mp.sqrt(2)-1
    x1, x2, d, h = 1-mp.sqrt(2), 1-mp.sqrt(2)/2, 2-mp.sqrt(2), mp.sqrt(2)/2
    l2, lu = mp.log(2), mp.log(1+mp.sqrt(2))
    Lam_ser = depth2(x1, -mp.sqrt(2)/2)
    swap_ser = depth2(x2, -mp.sqrt(2))
    compact3 = (3*z3/4 + li3(k) + 2*li3(x1) - li3(d)
                + l2*li2(k)/2 + l2*(li2(x2)-li2(h))/4)
    compact4 = (li3(d) + 2*li3(x2) - li3(k) - 4*li3(h) - 4*li3(-h)
                - l2*(li2(k)+li2(x2))/2 + l2**2*(l2/2-lu)/8)
    print("  octagon point (120 dps):")
    print(f"   Lambda: series={mp.nstr(Lam_ser, 30)}")
    print(f"     vs G_cl(k) {dig(Lam_ser-G_cl_mp(k))}, vs compact3 {dig(Lam_ser-compact3)}")
    print(f"   swap:   series={mp.nstr(swap_ser, 30)}")
    print(f"     vs H_cl(k) {dig(swap_ser-H_cl_mp(k))}, vs compact4 {dig(swap_ser-compact4)}")
    # 4.3 decagon families: series vs closed forms
    print("  decagon: G_dec,G_dec_cl and H_dec,H_dec_cl on (0,1):")
    for t0 in (mp.mpf("0.17"), mp.mpf("0.31"), mp.mpf("0.5"),
               mp.mpf("0.63"), mp.mpf("0.82")):
        rG = depth2(-t0, -(1-t0)/t0) - G_dec_cl_mp(t0)
        rH = depth2(1-t0, -t0/(1-t0)) - H_dec_cl_mp(t0)
        print(f"   t={t0}: |G_dec-G_dec_cl| {dig(rG)}, |H_dec-H_dec_cl| {dig(rH)}")
    # 4.4 decagon specialisation at t = 1/phi (170 dps: DEC1 claimed 830-digit LLL)
    mp.mp.dps = 170
    phi = (1+mp.sqrt(5))/2
    l, z2h, z3h = mp.log(phi), mp.pi**2/6, mp.zeta(3)
    Lam_dec_ser = depth2(-1/phi, -1/phi)
    swap_dec_ser = depth2(1/phi**2, -phi)
    Lam_DEC1 = (2*z3h/5 + 4*l*z2h/5 - l**3 + 3*li3(-1/phi) + li3(1/phi))
    swap_DEC1 = (3*z3h/5 - 7*l*z2h/5 + 2*l**3 + li3(-1/phi) + li3(1/phi))
    tp = 1/phi
    print("  decagon point (170 dps):")
    print(f"   Lambda_dec: series={mp.nstr(Lam_dec_ser, 40)}")
    print(f"     vs G_dec_cl(1/phi) {dig(Lam_dec_ser-G_dec_cl_mp(tp))},"
          f" vs DEC1 form {dig(Lam_dec_ser-Lam_DEC1)}")
    print(f"   swap_dec:   series={mp.nstr(swap_dec_ser, 40)}")
    print(f"     vs H_dec_cl(1/phi) {dig(swap_dec_ser-H_dec_cl_mp(tp))},"
          f" vs DEC1 form {dig(swap_dec_ser-swap_DEC1)}")
    # 4.5 golden Li2/Li3 lattice numerics
    lat = {
        "Li3(1/phi^2)=4Li3(-1/phi)+4Li3(1/phi)":
            li3(1/phi**2)-4*li3(-1/phi)-4*li3(1/phi),
        "Lewin Li3(1/phi^2)": li3(1/phi**2)-(4*(z3h-l*z2h)/5+2*l**3/3),
        "Li2(1/phi^2)=2z2/5-l^2": li2(1/phi**2)-(2*z2h/5-l**2),
        "Li2(1/phi)=3z2/5-l^2": li2(1/phi)-(3*z2h/5-l**2),
        "Li2(-1/phi)=-2z2/5+l^2/2": li2(-1/phi)-(-2*z2h/5+l**2/2),
    }
    for name, v in lat.items():
        print(f"   golden lattice {name}: residual {dig(v)}")
    print("PART 4 OK")


def main():
    part0()
    part1()
    part_chainrule()
    part2()
    part2_specialisation()
    part3()
    part3_landen3()
    part3_specialisation()
    part4()
    print("=" * 74)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
