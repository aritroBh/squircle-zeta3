#!/usr/bin/env python3
"""verify_all.py -- one command that re-checks the results of "Investigating zeta(3)".

This is an independent numerical re-check: it re-derives every quantity from its own integral
representation and compares it with the closed forms printed in the paper. It imports nothing from
the other folders of this repository. It only READS stored data for comparison:
  - squircle_square_end/outputs/M*_method{A,B}_1200.txt  (1200-digit M_0..M_3, methods A and B),
  - squircle_square_end/outputs/mk_deg9.tsv               (M_k, k <= 400, from an independent method).
Section "cert" RE-RUNS weight4_M3/w4_certified.py (python-flint ball arithmetic); that is a re-run of
the certified computation, not an independent check.

Run everything:
    python verify_all.py
or a subset of the sections:
    python verify_all.py squircle vform lemmas rhombus borel mk area corner hi cert pip ellipsoid cornerthm hexagon

Each check prints:  tag | residual | agreeing digits | digits required | PASS/FAIL, and a summary line
at the end. Digits = floor(-log10(|x-ref|/|ref|)) (relative), or -log10|x-ref| for identities "= 0".
"""
import sys
import time
from pathlib import Path
import mpmath as mp
import sympy as sp

HERE = Path(__file__).resolve().parent
T0 = time.time()
LOG = []


def digits(x, ref, mode):
    d = abs(x - ref)
    if d == 0:
        return 10**6
    if mode == "rel" and ref != 0:
        d = d / abs(ref)
    return int(mp.floor(-mp.log10(d)))


def check(tag, x, ref, need, mode="rel", note=""):
    dg = digits(x, ref, mode)
    ok = dg >= need
    res = x - ref
    dgs = ">=work" if dg == 10**6 else str(dg)
    line = f"  {tag:<58s} resid {mp.nstr(res, 3):>11s}  digits {dgs:>6s}  need {need:>4d}  {'PASS' if ok else '**FAIL**'}"
    if note:
        line += "  " + note
    print(line, flush=True)
    LOG.append((tag, dg, need, ok))
    return ok


def show(tag, x, n=40):
    print(f"  {tag:<40s} = {mp.nstr(x, n)}", flush=True)


def consts():
    r2 = mp.sqrt(2)
    C = dict(r2=r2, a=1 + r2, b=r2 - 1, L=mp.log(1 + r2), lam=mp.log(2), pi2=mp.pi**2,
             z3=mp.zeta(3))
    b = C["b"]
    C["Li2b"] = mp.polylog(2, b)
    C["Li3b"] = mp.polylog(3, b)
    C["Li3mb"] = mp.polylog(3, -b)
    C["Li3h"] = mp.polylog(3, 1 / r2)
    return C


# -------------------------------------------------------------------------------------------
# closed forms transcribed from the paper (Table I and Theorems 1-4)
# -------------------------------------------------------------------------------------------
def M0_closed(C):
    return C["lam"] - C["r2"] * C["L"]


def M1_closed(C):
    r2, L, lam, pi2 = C["r2"], C["L"], C["lam"], C["pi2"]
    return ((19 * r2 / 48 - mp.mpf(1) / 24) * pi2 + lam**2 / 2 + 2 * r2 * lam * L
            - (1 + 5 * r2 / 2) * L**2 - 5 * r2 * C["Li2b"])


def M1_closed_48(C):
    r2, L, lam, pi2 = C["r2"], C["L"], C["lam"], C["pi2"]
    return (-2 * pi2 + 19 * r2 * pi2 + 24 * lam**2 - 48 * L**2 - 120 * r2 * L**2
            + 96 * r2 * lam * L - 240 * r2 * C["Li2b"]) / 48


def T_closed(C):
    return (C["pi2"] + 6 * C["lam"]**2) / 24


def M2_closed(C):
    r2, L, lam, pi2, z3 = C["r2"], C["L"], C["lam"], C["pi2"], C["z3"]
    return (T_closed(C) + (189 * r2 / 64 - mp.mpf(19) / 16) * z3
            + (11 + 15 * r2) * C["Li3b"] + (11 - 15 * r2) * C["Li3mb"] - 27 * r2 * C["Li3h"]
            + (12 * L - 7 * r2 / 2 * lam) * C["Li2b"]
            + (15 * r2 / 8 - mp.mpf(7) / 12) * pi2 * L - (mp.mpf(1) / 24 + 65 * r2 / 48) * pi2 * lam
            + (mp.mpf(1) / 6 + 9 * r2 / 16) * lam**3 + 11 * r2 / 8 * lam**2 * L
            - (1 + 7 * r2 / 4) * lam * L**2 + (mp.mpf(25) / 6 - 7 * r2 / 3) * L**3)


def W3p_closed(C):
    r2, L, lam, pi2, z3 = C["r2"], C["L"], C["lam"], C["pi2"], C["z3"]
    return ((7 * r2 / 32 - mp.mpf(3) / 16) * z3 + (1 + r2) * C["Li3b"] + (1 - r2) * C["Li3mb"]
            - 2 * r2 * C["Li3h"] + (2 * L - r2 * lam) * C["Li2b"]
            + (r2 / 12 - mp.mpf(1) / 8) * pi2 * L - (1 + r2) / 24 * pi2 * lam + (1 + r2) / 24 * lam**3
            + r2 / 4 * lam**2 * L - (1 + r2) / 2 * lam * L**2 + (5 - r2) / 6 * L**3)


def conj4_terms(C):
    """Conjecture 4 as a list of (rational part, sqrt2 part, value); R = sum (p + q sqrt2) * value."""
    b, r2, L, lam, pi2, z3 = C["b"], C["r2"], C["L"], C["lam"], C["pi2"], C["z3"]
    Q = mp.mpf
    Li4 = lambda x: mp.polylog(4, x)
    Li3b, Li3mb, Li3h, Li2b = C["Li3b"], C["Li3mb"], C["Li3h"], C["Li2b"]
    t = [
        (1, -573, Li4(b)), (-1, -573, Li4(-b)), (0, 56, Li4(-b * b)),
        (0, 17, Li4(1 / r2)), (Q(-37) / 4, Q(-17) / 16, Li4(Q(1) / 2)),
        (10, 86, Li4(1 - 1 / r2)), (-10, 86, Li4(2 - r2)),
        (40, 7, L * Li3b), (-40, 7, L * Li3mb), (-72, 0, L * Li3h), (Q(63) / 8, Q(-3) / 16, L * z3),
        (6, 13, lam * Li3b), (6, -13, lam * Li3mb), (0, Q(-47) / 2, lam * Li3h),
        (Q(-3) / 8, Q(329) / 128, lam * z3),
        (30, 0, Li2b**2), (30, Q(7) / 2, Li2b * L**2), (-24, 0, Li2b * lam * L),
        (0, Q(-19) / 8, Li2b * lam**2), (Q(-19) / 4, Q(1) / 24, Li2b * pi2),
        (Q(103) / 384, Q(-4963) / 5760, pi2**2),
        (Q(67) / 24, Q(185) / 96, pi2 * L**2), (Q(-35) / 12, Q(41) / 12, pi2 * lam * L),
        (Q(5) / 12, Q(-1105) / 384, pi2 * lam**2),
        (Q(13) / 24, Q(-25) / 24, L**4), (Q(-73) / 6, Q(-29) / 3, L**3 * lam),
        (Q(63) / 8, Q(67) / 16, L**2 * lam**2), (Q(41) / 24, Q(17) / 6, L * lam**3),
        (Q(-61) / 192, Q(137) / 192, lam**4),
    ]
    return t


def conj4(C):
    t = conj4_terms(C)
    R0 = mp.fsum(p * v for p, q, v in t)
    R1 = mp.fsum(q * v for p, q, v in t)
    return R0 + C["r2"] * R1, R0, R1


# -------------------------------------------------------------------------------------------
# Section 1: squircle coefficients M_0..M_4 from the t-form, by symbolic differentiation in eps
# -------------------------------------------------------------------------------------------
def squircle_Mk(K):
    """M_k = -(1/k!) int_0^{1/2} d^k/deps^k Phi_eps(t) |_{eps=0} dt, with the cancellation-free form
    Phi = b (1 - z/(1+sqrt(1+z^2))), b = (1-t)^(eps-1), z = (t/(1-t))^(1-eps)  (t < 1/2)."""
    e, T1, T2 = sp.symbols("e T1 T2")
    bb = sp.exp((e - 1) * T1)
    z = sp.exp((1 - e) * T2)
    expr = bb * (1 - z / (1 + sp.sqrt(1 + z**2)))
    fns = []
    for k in range(K + 1):
        fns.append(sp.lambdify((T1, T2), expr.subs(e, 0), "mpmath"))
        expr = sp.diff(expr, e)
    out = []
    for k in range(K + 1):
        f = fns[k]
        g = lambda t, f=f: f(mp.log1p(-t), mp.log(t) - mp.log1p(-t))
        val = mp.quad(g, [0, mp.mpf(1) / 1024, mp.mpf(1) / 64, mp.mpf(1) / 8, mp.mpf(1) / 2])
        out.append(-val / mp.factorial(k))
    return out


M3_NOTE = "0.02277052149178657745673050974232780025018369172921726137909"   # note, 58 digits (string: convert at use)
W4_NOTE = "0.795785059832240811320359819629726762546527413871"               # note, 48 digits
M2_NOTE = "-0.022161911087269399423754937739791514977955"
KAPPA_NOTE = "2.02095276965635835694690816772523888426901720579"
M4_STORED = "0.0060279931126091225045161082508776055714767341482928028954590935391901689357716143"


def sec_squircle():
    print("\n== 1. Squircle coefficients M_0..M_4 (t-form, symbolic eps-derivatives, dps 70) ==")
    mp.mp.dps = 70
    C = consts()
    t = time.time()
    M = squircle_Mk(4)
    print(f"  [computed M_0..M_4 in {time.time() - t:.0f}s]")
    for k in range(5):
        show(f"M_{k}", M[k], 60)
    check("M0  M0 = log2 - sqrt2 L", M[0], M0_closed(C), 60)
    check("Thm1  M1 = Theorem 1 (expanded form)", M[1], M1_closed(C), 60)
    check("Thm1  M1 = Theorem 1 (48 M1 form)", M[1], M1_closed_48(C), 60)
    check("Thm1  kappa = 8 M1 vs printed value", 8 * M[1], mp.mpf(KAPPA_NOTE), 45)
    check("Thm1  Ahmed: kappa = 2.0209527696...", 8 * M[1], mp.mpf("2.0209527696"), 10)
    check("Thm2  M2 = Theorem 2", M[2], M2_closed(C), 60)
    check("Thm2  M2 vs printed value", M[2], mp.mpf(M2_NOTE), 40)
    check("Thm3  M3 vs printed 58-digit value", M[3], mp.mpf(M3_NOTE), 55)
    check("Prop1  M4 vs stored table (independent method)", M[4], mp.mpf(M4_STORED), 55)
    return M


# -------------------------------------------------------------------------------------------
# Section 2: v-form integrals T, W3, W3', W4, U, V and the M3 split, Conjecture 4
# -------------------------------------------------------------------------------------------
def sec_vform(M=None):
    print("\n== 2. v-form integrals, the M3 split and the closed form R of W4 (dps 70) ==")
    mp.mp.dps = 70
    C = consts()
    a, b, r2 = C["a"], C["b"], C["r2"]
    if M is None:
        M = squircle_Mk(3)
    ell = lambda v: mp.log(v - b) + mp.log(v + a) - mp.log(4 * v)
    em = lambda v: mp.log(a - v) + mp.log(v + b) - mp.log(4 * v)
    vq = lambda f: mp.quad(f, [b, mp.mpf(1) / 2, 1, 2, a])
    M1v = vq(lambda v: ell(v) * ((1 - 1 / r2) / v - 2 / (v + a)))
    check("Thm1  M1 from eq.(M1v) quadrature", M1v, M1_closed(C), 60)
    T = vq(lambda v: ell(v) * (1 / (2 * v) + 1 / (2 * (v - a)) - 1 / (v + a) - 1 / (2 * (v + b))))
    check("Thm2  T = (pi^2 + 6 log^2 2)/24  [Prop. T]", T, T_closed(C), 60)
    check("Thm2  T = (zeta(2) - Li2(1/2))/2", T_closed(C), (mp.zeta(2) - mp.polylog(2, mp.mpf(1) / 2)) / 2, 60)
    W3 = vq(lambda v: ell(v)**2 * ((4 - 3 * r2) / (8 * v) - 1 / (v + a)) + r2 / 8 * ell(v) * em(v) / v)
    check("Thm2  M2 = W3 + T (quadrature of eq.(W3))", W3 + T, M[2], 60)
    W3p = vq(lambda v: ell(v)**2 * ((1 - r2) / (4 * v) + 1 / (4 * (v - a)) - 1 / (4 * (v + a)) - 1 / (2 * (v + b)))
             + ell(v) * em(v) * ((1 + r2) / (4 * v) - 1 / (2 * (v + a))))
    check("Thm3  W3' closed form vs quadrature", W3p, W3p_closed(C), 60)
    U = vq(lambda v: ell(v)**3 * (1 / (6 * v) - 1 / (3 * (v + a))))
    V = vq(lambda v: -mp.mpf(7) / 48 * ell(v)**3 / v + ell(v)**2 * em(v) / (16 * v))
    W4 = U + r2 * V
    W4b = vq(lambda v: ell(v)**3 * ((8 - 7 * r2) / (48 * v) - 1 / (3 * (v + a))) + r2 / 16 * ell(v)**2 * em(v) / v)
    check("Thm3  W4 = U + sqrt2 V (two quadratures agree)", W4, W4b, 60)
    check("Thm3  M3 = W4 + W3' + T  [Theorem 3 split]", W4 + W3p_closed(C) + T_closed(C), M[3], 60)
    check("Thm3  W4 vs printed value", W4, mp.mpf(W4_NOTE), 45)
    R, R0, R1 = conj4(C)
    check("Thm4  W4 = Conjecture 4 RHS", W4, R, 60)
    check("Thm4  U - R0 = 0  (rational-prefactor part)", U - R0, mp.mpf(0), 60, "abs")
    check("Thm4  V - R1 = 0  (sqrt2-prefactor part)", V - R1, mp.mpf(0), 60, "abs")
    show("(W4 - R)/pi^4", (W4 - R) / mp.pi**4, 5)


# -------------------------------------------------------------------------------------------
# Section 3: dilog / trilog lemmas, eq.(chi3), other forms of M2
# -------------------------------------------------------------------------------------------
def sec_lemmas():
    print("\n== 3. Lemma 2 (Li2), the Li3 identities of Appendix B, and two other forms of M2 (dps 80) ==")
    mp.mp.dps = 80
    C = consts()
    b, r2, L, lam, pi2, z3 = C["b"], C["r2"], C["L"], C["lam"], C["pi2"], C["z3"]
    Li2, Li3 = (lambda x: mp.polylog(2, x)), (lambda x: mp.polylog(3, x))
    Li2b, Li3b, Li3mb, Li3h = C["Li2b"], C["Li3b"], C["Li3mb"], C["Li3h"]
    dil = [
        ("Li2(b)-Li2(-b)", Li2(b) - Li2(-b), pi2 / 8 - L**2 / 2),
        ("Li2(b^2)", Li2(b * b), 4 * Li2b - pi2 / 4 + L**2),
        ("Li2(2-sqrt2)", Li2(2 - r2), -Li2b + pi2 / 6 + lam * L / 2 - L**2),
        ("Li2(-sqrt2)", Li2(-r2), Li2b - pi2 / 6 - lam * L / 2 + L**2 / 2),
        ("Li2(-b^2)", Li2(-b * b), 6 * Li2b - 11 * pi2 / 24 + 2 * L**2),
        ("Li2((2+sqrt2)/4)", Li2((2 + r2) / 4), 6 * Li2b - 7 * pi2 / 24 - 9 * lam**2 / 8 - 3 * lam * L / 2 + 7 * L**2 / 2),
    ]
    for n, x, y in dil:
        check("Lemma2 " + n, x, y, 70)
    Q = mp.mpf
    tri = [
        ("Li3(-a)", Li3(-C["a"]), Li3mb - pi2 * L / 6 - L**3 / 6),
        ("Li3(1-1/sqrt2)", Li3(1 - 1 / r2), z3 - Li3mb - Li3h - pi2 * lam / 12 + lam**3 / 24 + lam**2 * L / 8),
        ("Li3(-1/sqrt2)", Li3(-1 / r2), Q(7) / 32 * z3 - Li3h - pi2 * lam / 48 + lam**3 / 24),
        ("Li3(2-sqrt2)", Li3(2 - r2), Q(25) / 32 * z3 - Li3b + Li3h + 5 * pi2 * lam / 48 - pi2 * L / 6
         - lam**3 / 48 - lam * L**2 / 4 + L**3 / 3),
        ("Li3(-b^2)", Li3(-b * b), -Q(11) / 8 * z3 + 6 * Li3b + 6 * Li3mb + pi2 * L / 6 - L**3 / 3),
        ("Li3(2b)", Li3(2 * b), -Q(3) / 8 * z3 - 6 * Li3b + 2 * Li3mb + 8 * Li3h + 7 * pi2 * lam / 12
         - 2 * pi2 * L / 3 - lam**3 / 6 - lam * L**2 + 4 * L**3 / 3),
        ("Li3(-b/2)", Li3(-b / 2), Q(11) / 8 * z3 + 2 * Li3b - 6 * Li3mb - 8 * Li3h - 5 * pi2 * lam / 12
         + pi2 * L / 2 + lam**3 / 3 + lam**2 * L / 2 - lam * L**2 / 2 - L**3 / 2),
        ("Li3(-2a)", Li3(-2 * C["a"]), Q(11) / 8 * z3 + 2 * Li3b - 6 * Li3mb - 8 * Li3h - 7 * pi2 * lam / 12
         + pi2 * L / 3 + lam**3 / 6 - lam * L**2 - 2 * L**3 / 3),
        ("Li3((2-sqrt2)/4)", Li3((2 - r2) / 4), Q(101) / 32 * z3 + 7 * Li3b - 13 * Li3mb - 18 * Li3h
         - 17 * pi2 * lam / 16 + 7 * pi2 * L / 6 + 15 * lam**3 / 16 + 9 * lam**2 * L / 8 - 3 * lam * L**2 / 4
         - 4 * L**3 / 3),
        ("Li3((2+sqrt2)/4)", Li3((2 + r2) / 4), -Q(25) / 32 * z3 - 13 * Li3b + 7 * Li3mb + 18 * Li3h
         + 13 * pi2 * lam / 16 - 7 * pi2 * L / 6 + 3 * lam**3 / 16 - 9 * lam**2 * L / 8 - 3 * lam * L**2 / 4
         + 7 * L**3 / 3),
    ]
    for n, x, y in tri:
        check("AppB-Li3 " + n, x, y, 70)
    chi8 = {1: 1, 3: -1, 5: -1, 7: 1}
    L8 = mp.fsum(s * mp.zeta(3, mp.mpf(r) / 8) for r, s in chi8.items()) / 512
    chi3b = (Li3b - Li3mb) / 2
    check("M2-forms eq.(chi3): 48 chi3(b) = 32sqrt2 L(chi8,3) - 3pi^2 L + 4L^3",
          48 * chi3b, 32 * r2 * L8 - 3 * pi2 * L + 4 * L**3, 70)
    Cl3 = mp.re(mp.polylog(3, mp.expjpi(mp.mpf(1) / 4)))
    check("M2-forms Cl3(pi/4) = (sqrt2/2) L(chi8,3) - 3 zeta(3)/256", Cl3, r2 / 2 * L8 - 3 * z3 / 256, 70)
    M2 = M2_closed(C)
    l2 = lam
    formB1 = T_closed(C) + (-228 * z3 + 4224 * Li3b + 7680 * L8 - 8 * pi2 * l2 + 152 * pi2 * L + 32 * l2**3
                            - 192 * l2 * L**2 + 448 * L**3 + 2304 * L * Li2b
                            + r2 * (567 * z3 - 5184 * Li3h - 2816 * L8 - 260 * pi2 * l2 + 108 * l2**3
                                    + 264 * l2**2 * L - 336 * l2 * L**2 + 32 * L**3 - 672 * l2 * Li2b)) / 192
    check("M2-forms App.B form with L(chi8,3) = Theorem 2", formB1, M2, 70)
    formB2 = (pi2 / 24 + l2**2 / 4 + l2**3 / 6 - l2 * L**2 + Q(7) / 3 * L**3 - pi2 * l2 / 24 + Q(19) / 24 * pi2 * L
              + 12 * L * Li2b - Q(49) / 32 * z3 - Q(88) / 3 * Cl3 + 22 * Li3b
              + r2 * (Q(9) / 16 * l2**3 + Q(11) / 8 * l2**2 * L - Q(7) / 4 * l2 * L**2 + L**3 / 6
                      - Q(65) / 48 * pi2 * l2 - Q(7) / 2 * l2 * Li2b + Q(219) / 64 * z3 + 40 * Cl3 - 27 * Li3h))
    check("M2-forms App.B form with Cl3(pi/4) = Theorem 2", formB2, M2, 70)
    W4rem = W3p_closed(C) - z3 / 8 - pi2 * lam / 24 - lam**3 / 12
    show("Prop1 weight-3 part of M4 (Remark 12) =", W4rem, 30)


# -------------------------------------------------------------------------------------------
# Section 4: rhombus end, general semi-axes
# -------------------------------------------------------------------------------------------
def P_direct(p, a, b):
    """Perimeter of |x/a|^p + |y/b|^p = 1 from x = a t^e, y = b (1-t)^e, e = 1/p (no other input)."""
    e = 1 / p
    a, b = mp.mpf(a), mp.mpf(b)
    f = lambda t: mp.sqrt(a**2 * t**(2 * e - 2) + b**2 * (1 - t)**(2 * e - 2))
    return 4 * e * mp.quad(f, [0, mp.mpf(1) / 64, mp.mpf(1) / 2, mp.mpf(63) / 64, 1])


def P_R2rep(p, a, b):
    """R2's exact representation P = 4r(1+eta) int [t(1-t)]^(eta/2) g(eta l) dt."""
    a, b = mp.mpf(a), mp.mpf(b)
    r = mp.sqrt(a * a + b * b)
    q = (a * a - b * b) / (a * a + b * b)
    eta = 1 / p - 1
    g = lambda z: mp.sqrt(mp.cosh(z) + q * mp.sinh(z))
    f = lambda t: (t * (1 - t))**(eta / 2) * g(eta * mp.log(t / (1 - t)))
    return 4 * r * (1 + eta) * mp.quad(f, [0, mp.mpf(1) / 64, mp.mpf(1) / 2, mp.mpf(63) / 64, 1])


def sec_rhombus():
    print("\n== 4. Rhombus end p -> 1, general semi-axes (dps 95, 9-point stencil, h = 1e-8) ==")
    mp.mp.dps = 95
    z3, pi2 = mp.zeta(3), mp.pi**2
    offs = list(range(-4, 5))
    Wt = sp.finite_diff_weights(4, offs, 0)
    W = [[mp.mpf(sp.Rational(w).p) / sp.Rational(w).q for w in Wt[m][-1]] for m in range(5)]
    h = mp.mpf(10)**-8
    for (a, b) in [(1, 1), (1, 2), (3, 1), (2, 5), (1, 7)]:
        a, b = mp.mpf(a), mp.mpf(b)
        r = mp.sqrt(a * a + b * b)
        t = time.time()
        vals = [P_direct(1 + j * h, a, b) for j in offs]
        D = [mp.fsum(w * v for w, v in zip(W[m], vals)) / h**m for m in range(5)]
        if a == b:
            # Prop. 2: P/sqrt2 = sum a_n d^n, a_4 = -3 zeta3 - 11 pi^4/480 + pi^2/2  (4th derivative, same stencil)
            check("Prop2 rhombus a_4 from direct 4th derivative (square)", D[4] / (24 * mp.sqrt(2)),
                  -3 * z3 - 11 * mp.pi**4 / 480 + pi2 / 2, 25)
        ex = [4 * r, mp.mpf(0), 4 * pi2 * a * a * b * b / (3 * r**3), 24 * a * a * b * b / r**3 * (z3 - pi2 / 3)]
        tag = f"({int(a)},{int(b)})"
        check(f"Thm5 P(1)={tag}: 4r", vals[4], ex[0], 80)
        check(f"Thm5 P'(1)={tag}: 0", D[1], ex[1], 45, "abs")
        check(f"Thm5 P''(1)={tag}: 4pi^2a^2b^2/(3r^3)", D[2], ex[2], 45)
        check(f"Thm5 P'''(1)={tag}: 24a^2b^2/r^3 (zeta3 - pi^2/3)", D[3], ex[3], 35,
              note=f"[{time.time() - t:.0f}s]")
        if a == b:
            r2 = mp.sqrt(2)
            check("Thm5 square: P''(1) = (sqrt2/3) pi^2", D[2], r2 / 3 * pi2, 45)
            check("Thm5 square: P'''(1) = 6sqrt2 zeta3 - 2sqrt2 pi^2", D[3], 6 * r2 * z3 - 2 * r2 * pi2, 35)
            check("RM eq.(104) coefficient: P''(1)/2 = 4r pi^2 a^2b^2/(6 r^4)", D[2] / 2,
                  4 * r * pi2 * a * a * b * b / (6 * r**4), 45)
    mp.mp.dps = 50
    for p in [mp.mpf("0.6"), mp.mpf("1.3"), mp.mpf(3)]:
        check(f"Eq15 exact representation vs direct, (a,b)=(1,2), p={p}", P_R2rep(p, 1, 2), P_direct(p, 1, 2), 45)
    hh = lambda t: mp.log(t * (1 - t)) / 2
    ll = lambda t: mp.log(t / (1 - t))
    q = lambda f: mp.quad(f, [0, mp.mpf(1) / 2, 1])
    check("Eq15 moment int h^3 = -6 + 3zeta3/2 + pi^2/4", q(lambda t: hh(t)**3), -6 + 3 * mp.zeta(3) / 2 + mp.pi**2 / 4, 45)
    check("Eq15 moment int h l^2 = -(2zeta3 + pi^2/3)", q(lambda t: hh(t) * ll(t)**2), -(2 * mp.zeta(3) + mp.pi**2 / 3), 45)
    check("Eq15 moment int h^2 = 2 - pi^2/12", q(lambda t: hh(t)**2), 2 - mp.pi**2 / 12, 45)


# -------------------------------------------------------------------------------------------
# Section 5: zero radius at p = 1 -- growth of the formal Taylor coefficients (square)
# -------------------------------------------------------------------------------------------
def sec_borel(N=64):
    print(f"\n== 5. Zero radius at the rhombus: Taylor coefficients S_n of P/(4r) in eta, n <= {N} (dps 90) ==")
    mp.mp.dps = 90
    q = mp.mpf(0)   # square
    # c_k: Taylor coefficients of g(z) = sqrt(cosh z + q sinh z)
    ser = [(1 if k % 2 == 0 else q) / mp.factorial(k) for k in range(N + 1)]
    c = [mp.mpf(1)]
    for n in range(1, N + 1):
        c.append((ser[n] - mp.fsum(c[k] * c[n - k] for k in range(1, n))) / 2)
    # homogeneous parts of log A(eta, s); index i = power of s
    g = [[mp.mpf(0)] * (m + 1) for m in range(N + 1)]
    for m in range(2, N + 1):
        zm = (-1)**m * mp.zeta(m) / m
        for i in range(m + 1):
            v = mp.binomial(m, i) * mp.mpf(2)**(i - m) * (1 + (-1)**i)
            if i == 0:
                v -= 1
            g[m][i] = zm * v
    f = [[mp.mpf(1)]]
    for n in range(1, N + 1):
        acc = [mp.mpf(0)] * (n + 1)
        for j in range(2, n + 1):
            fj = f[n - j]
            for i1, x in enumerate(g[j]):
                if x == 0:
                    continue
                for i2, y in enumerate(fj):
                    acc[i1 + i2] += j * x * y
        f.append([x / n for x in acc])
    S = [mp.fsum(c[k] * mp.factorial(k) * f[n][k] for k in range(0, n + 1, 2)) for n in range(N + 1)]
    check("Thm6 S_2 = pi^2/24 (square)", S[2], mp.pi**2 / 24, 50)
    check("Thm6 S_3 = -zeta(3)/4 (square)", S[3], -mp.zeta(3) / 4, 50)
    # rho_n = -S_{n+2} / ((n+1)(n+2) S_n) -> 1/R^2 ; Richardson in 1/n over each parity
    print("   n   -S_{n+2}/((n+1)(n+2)S_n)      R estimate (Richardson order 3)")
    for par in (0, 1):
        ns = [n for n in range(20 + par, N - 1, 2)]
        rho = {n: -S[n + 2] / ((n + 1) * (n + 2) * S[n]) for n in ns}
        for n in ns[-4:]:
            # Richardson on rho_n = x + a/n + b/n^2 + c/n^3 using n, n-2, n-4, n-6
            pts = [n - 6, n - 4, n - 2, n]
            A = mp.matrix([[1, mp.mpf(1) / k, mp.mpf(1) / k**2, mp.mpf(1) / k**3] for k in pts])
            y = mp.matrix([rho[k] for k in pts])
            x = mp.lu_solve(A, y)[0]
            print(f"  {n:3d}  {mp.nstr(rho[n], 15):>22s}   R = {mp.nstr(1 / mp.sqrt(x), 12)}")
        Rest = 1 / mp.sqrt(x)
        check(f"Thm6 Borel radius R = pi/2 from S_n growth (parity {par})", Rest, mp.pi / 2, 4)
    show("|S_n/n!|^(-1/n) at n = N", abs(S[N] / mp.factorial(N))**(-mp.mpf(1) / N), 12)

    # ---- Prop. 2: a_n of P/sqrt2 = sum a_n d^n (d = p-1), from S via eta = -d/(1+d):
    #      [d^n] eta^m = (-1)^n C(n-1, m-1)  (n >= m >= 1)
    a = [4 * S[0]] + [4 * mp.fsum((-1)**n * mp.binomial(n - 1, m - 1) * S[m] for m in range(1, n + 1))
                      for n in range(1, N + 1)]
    z3, z5, z7, pi = mp.zeta(3), mp.zeta(5), mp.zeta(7), mp.pi
    Q = mp.mpf
    d11 = {  # transcribed from Proposition 2 and Table I of the paper
        2: pi**2 / 6,
        3: z3 - pi**2 / 3,
        4: -3 * z3 - 11 * pi**4 / 480 + pi**2 / 2,
        5: -7 * pi**2 * z3 / 24 + 6 * z3 - 3 * z5 / 4 - 2 * pi**2 / 3 + 11 * pi**4 / 120,
        6: -7 * z3**2 / 8 - 10 * z3 + 35 * pi**2 * z3 / 24 + 15 * z5 / 4 - 11 * pi**4 / 48 + 5 * pi**2 / 6
           + 577 * pi**6 / 34560,
        7: 21 * z3**2 / 4 - 35 * pi**2 * z3 / 8 + 15 * z3 + 119 * pi**4 * z3 / 640 - 45 * z5 / 4
           + 47 * pi**2 * z5 / 32 + 129 * z7 / 16 - 577 * pi**6 / 5760 - pi**2 + 11 * pi**4 / 24,
        8: -147 * z3**2 / 8 + 169 * pi**2 * z3**2 / 192 + 141 * z3 * z5 / 16 - 833 * pi**4 * z3 / 640 - 21 * z3
           + 245 * pi**2 * z3 / 24 - 329 * pi**2 * z5 / 32 + 105 * z5 / 4 - 903 * z7 / 16
           - Q(2891923) * pi**8 / 116121600 - 77 * pi**4 / 96 + 7 * pi**2 / 6 + 4039 * pi**6 / 11520,
    }
    check("Prop2 rhombus a_0 = 4", a[0], mp.mpf(4), 80)
    check("Prop2 rhombus a_1 = 0", a[1], mp.mpf(0), 80, "abs")
    for n, v in d11.items():
        check(f"Prop2 rhombus a_{n} closed form vs own Gamma generating function", a[n], v, 80)
    check("Prop2 rhombus printed a_8 = -160.30725951244082378376571501688937243223996019071",
          a[8], mp.mpf("-160.30725951244082378376571501688937243223996019071"), 48)
    # f_n = Taylor coefficients of F = P/(4 sqrt2 (1+eta)) = S/(1+eta)
    f = [mp.fsum((-1)**(n - j) * S[j] for j in range(n + 1)) for n in range(N + 1)]
    check("Prop2 rhombus printed f_60 (65 digits)", f[60],
          mp.mpf("-3.0158042558052800566810935892480634882385393222591623494933537906e+67"), 60)
    # growth law f_n ~ sqrt2 cos(pi(3-2n)/4) n! (2/pi)^n n^(-3/2): Richardson in 1/n over n = n0, n0+4, ...
    Cn = {n: f[n] * mp.mpf(n)**1.5 / (mp.factorial(n) * (2 / pi)**n) / (mp.sqrt(2) * mp.cos(pi * (3 - 2 * n) / 4))
          for n in range(30, N + 1)}
    print("   n    C_n = f_n n^1.5/(n!(2/pi)^n sqrt2 cos(pi(3-2n)/4))   (expected: C_60 = 0.985479961182, C_59 = 1.02641818992)")
    for n in [56, 57, 58, 59, 60]:
        print(f"  {n:3d}  {mp.nstr(Cn[n], 14)}")
    check("Prop2 rhombus printed C_60 = 0.985479961182", Cn[60], mp.mpf("0.985479961182"), 11)
    check("Prop2 rhombus printed C_59 = 1.02641818992", Cn[59], mp.mpf("1.02641818992"), 10)
    for pts in ([48, 52, 56, 60], [49, 53, 57, 61], [50, 54, 58, 62], [51, 55, 59, 63]):
        A = mp.matrix([[1, mp.mpf(1) / k, mp.mpf(1) / k**2, mp.mpf(1) / k**3] for k in pts])
        K = mp.lu_solve(A, mp.matrix([Cn[k] for k in pts]))[0]
        check(f"Prop2 growth-law constant K -> 1 (Richardson, n = {pts})", K, mp.mpf(1), 4,
              note=f"K = {mp.nstr(K, 9)}")


# -------------------------------------------------------------------------------------------
# Section 6: high-order M_k by a Cauchy integral (u-form), vs asymptotic formulas
# -------------------------------------------------------------------------------------------
def M_of_q(q):
    """M(q) = -(1/4) int_0^inf (1+e^-u)^(1-q) [1 - X/(1+sqrt(1+X^2))] / cosh^2(u/2) du, X = e^{-(1-q)u};
    valid (holomorphic) for Re q < 1.  Derived independently from P/8 = 1 + q M(q), q = 1/p."""
    def f(u):
        X = mp.exp(-(1 - q) * u)
        return (1 + mp.exp(-u))**(1 - q) * (1 - X / (1 + mp.sqrt(1 + X * X))) / mp.cosh(u / 2)**2
    return -mp.quad(f, [0, 1, 4, 16, 40, 80, 160, mp.inf]) / 4


def sec_mk(N=512, rad="0.9", dps=50):
    # Aliasing error of the N-point rule is about M_{k+N} r^N; since |M_k| ~ e^{-sqrt(pi k)}/k decays slowly,
    # r = 0.97 leaves an error ~2e-27 (seen in a test run).  r = 0.9 makes aliasing ~1e-44 while the
    # amplification 1/r^k of quadrature error stays below 1e12 for k <= 250.
    print(f"\n== 6. M_k for large k: Cauchy integral on |q| = {rad}, N = {N} nodes (dps {dps}) ==")
    mp.mp.dps = dps
    C = consts()
    r = mp.mpf(rad)
    t = time.time()
    vals = {}
    for j in range(N // 2 + 1):
        vals[j] = M_of_q(r * mp.expjpi(mp.mpf(2 * j) / N))
    for j in range(N // 2 + 1, N):
        vals[j] = mp.conj(vals[N - j])
    print(f"  [{N // 2 + 1} complex quadratures in {time.time() - t:.0f}s]")

    def coef(k):
        s = mp.fsum(vals[j] * mp.expjpi(-mp.mpf(2 * j * k) / N) for j in range(N))
        return mp.re(s) / N / r**k
    check("Lem1  Cauchy M0 vs closed form", coef(0), M0_closed(C), 30)
    check("Lem1  Cauchy M1 vs Theorem 1", coef(1), M1_closed(C), 30)
    check("Lem1  Cauchy M2 vs Theorem 2", coef(2), M2_closed(C), 30)
    check("Lem1  Cauchy M3 vs printed value", coef(3), mp.mpf(M3_NOTE), 30)
    d13b = {}
    try:
        fn = HERE / "squircle_square_end" / "outputs" / "mk_deg9.tsv"
        for line in open(fn):
            k, v, _ = line.rstrip("\n").split("\t")
            d13b[int(k)] = mp.mpf(v.split(" ")[0].strip("["))
    except OSError:
        pass
    print("    k   M_k (Cauchy)            law -(sqrt(pi)/k)e^{-sqrt(pi k)}cos   ratio   (sqrt(pi)/(2k))..cos ratio")
    spi = mp.sqrt(mp.pi)
    for k in [10, 50, 100, 150, 200, 225, 250]:
        mk = coef(k)
        x = mp.sqrt(mp.pi * k)
        law = -(spi / k) * mp.exp(-x) * mp.cos(x)
        d13 = (spi / (2 * k)) * mp.exp(-x) * mp.cos(x)
        print(f"  {k:4d}  {mp.nstr(mk, 14):>22s}   {mp.nstr(law, 8):>14s}   {mp.nstr(mk / law, 6):>8s}"
              f"   {mp.nstr(mk / d13, 6):>8s}")
        if k in d13b:
            check(f"Obs1  M_{k}: Cauchy vs stored table", mk, d13b[k], 12 if k <= 200 else 9)
    # amplitude/phase test away from zeros of cos
    for k in [150, 200, 250]:
        mk = coef(k)
        x = mp.sqrt(mp.pi * k)
        check(f"Obs1  M_{k} / asymptotic law (16) -> 1 (O(k^-1/2))", mk, -(spi / k) * mp.exp(-x) * mp.cos(x), 1,
              note="leading-order only")


# -------------------------------------------------------------------------------------------
# Section 7: area of z = x^n y^n over [0,1]^2
# -------------------------------------------------------------------------------------------
def area_minus_1(n, deg=5):
    """A(n) - 1 with x = e^{-g/n}, y = e^{-(w-g)/n} (derived here):
    A - 1 = n^-2 int_0^inf e^{-w/n} int_0^w [sqrt(1 + n^2 e^{-2w}(e^{2g/n} + e^{2(w-g)/n})) - 1] dg dw.
    The inner integrand is analytic and slowly varying in g (scale n), so a fixed Gauss-Legendre rule
    (3*2^(deg-1) nodes) is used for it; tanh-sinh for the outer integral."""
    n = mp.mpf(n)
    nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(deg, mp.mp.prec)

    def inner(w):
        h = w / 4       # map [-1,1] -> [0, w/2]
        s = mp.fsum(wt * (mp.sqrt(1 + n * n * mp.exp(-2 * w) * (mp.exp(2 * (h * (x + 1)) / n)
                                                                + mp.exp(2 * (w - h * (x + 1)) / n))) - 1)
                    for x, wt in nodes)
        return 2 * h * s * mp.exp(-w / n)
    c = mp.log(n)
    pts = [0, max(c - 6, mp.mpf(1)) / 2, max(c - 6, mp.mpf(1)), c, c + 4, c + 12, 2 * c + 40, mp.inf]
    pts = sorted(set(pts))
    return mp.quad(inner, pts) / n**2


def area_direct(n):
    """Brute-force 2D quadrature in the original variables (for a small n only)."""
    f = lambda x, y: mp.sqrt(1 + n * n * x**(2 * n - 2) * y**(2 * n - 2) * (x * x + y * y)) - 1
    br = [0, mp.mpf("0.5"), mp.mpf("0.8"), mp.mpf("0.95"), 1]
    return mp.quad(f, br, br)


def sec_area():
    print("\n== 7. Area of z = x^n y^n over the unit square (dps 30) ==")
    mp.mp.dps = 30
    r2 = mp.sqrt(2)
    mp.mp.dps = 22
    check("Prop3 reduction formula vs brute-force 2D quadrature, n = 5", area_minus_1(5), area_direct(5), 15)
    mp.mp.dps = 30
    check("Prop3 inner Gauss-Legendre rule: 48 vs 96 nodes, n = 1000", area_minus_1(1000, 5), area_minus_1(1000, 6), 25)
    print("        n        A(n) - 1                 n(A-1) - sqrt2        R(n) = n^2(A-1-sqrt2/n)   R(n)+log^2(n)/2")
    data = {}
    for n in [10, 100, 1000, 10**4, 10**5, 10**6]:
        t = time.time()
        E = area_minus_1(n)
        Rn = n * n * (E - r2 / n)
        data[n] = Rn
        print(f"  {n:>8d}  {mp.nstr(E, 22):>24s}  {mp.nstr(n * E - r2, 8):>14s}  {mp.nstr(Rn, 12):>20s}"
              f"  {mp.nstr(Rn + mp.log(n)**2 / 2, 10):>14s}   [{time.time() - t:.0f}s]")
    check("Prop3 n^2(A-1-sqrt2/n) at n = 1e5 = -76.7411... (25-digit run)", data[10**5], mp.mpf("-76.7411"), 5)
    # elementary bounds (Proposition 3 of the paper):  lo(n) <= R(n) <= up(n) <= 0
    print("   n     lower bound          R(n)              upper bound   (bounds of Proposition 3)")
    for n in sorted(data):
        nn = mp.mpf(n)
        c = mp.log(r2 * nn)
        lo = -(c**2 / 2 + c + 1) * nn**2 / (nn - 1)**2 - 2 * r2 * nn / (nn + 1)
        g0 = mp.log(nn)
        cp = mp.log(r2 * nn / g0)
        up = -(1 - 1 / (2 * g0)) * mp.exp(-cp / (nn - mp.mpf(1) / 2)) * cp**2 / 2 * nn**2 / (nn - mp.mpf(1) / 2)**2
        ok = lo <= data[n] <= up <= 0
        print(f"  {n:>8d}  {mp.nstr(lo, 10):>14s}  {mp.nstr(data[n], 12):>16s}  {mp.nstr(up, 10):>14s}   "
              f"{'PASS' if ok else '**FAIL**'}")
        LOG.append((f"Prop3 bound lo <= R(n) <= up, n = {n}", 0 if ok else -1, 0, ok))
        print(f"            R(n)/log^2 n = {mp.nstr(data[n] / g0**2, 6)}   (limit -1/2 from the F14 bounds)")
    # beta, gamma fit from the 3 largest n: R + log^2/2 = beta log n + gamma + delta log^2(n)/n
    ns = [10**4, 10**5, 10**6]
    A = mp.matrix([[mp.log(n), 1, mp.log(n)**2 / n] for n in ns])
    y = mp.matrix([data[n] + mp.log(n)**2 / 2 for n in ns])
    sol = mp.lu_solve(A, y)
    print(f"  fit R(n) + log^2(n)/2 ~ beta log n + gamma (+ d log^2 n / n): beta = {mp.nstr(sol[0], 8)}, "
          f"gamma = {mp.nstr(sol[1], 8)}  (DATA ONLY; not a result)")
    show("n(A-1) - sqrt2 at n = 1e6", 10**6 * (data[10**6] / mp.mpf(10)**12 + r2 / 10**6) - r2, 6)


# -------------------------------------------------------------------------------------------
# Section 8: corner coefficients M_k(alpha) at interior angle alpha
# -------------------------------------------------------------------------------------------
def corner_coeffs(alpha, K):
    """M(u;alpha) = int_0^1 [(1+g)^(-1-u) sqrt(1 + g^(2-2u) - 2cos(alpha) g^(1-u)) - 1] g^(u-1) dg
    = sum_k M_k(alpha) u^k ;  M_k by symbolic u-derivatives at u = 0 (own implementation)."""
    u, g, Lg, Lp, c = sp.symbols("u g Lg Lp c")
    expr = (sp.exp(-(1 + u) * Lp) * sp.sqrt(1 + g**2 * sp.exp(-2 * u * Lg) - 2 * c * g * sp.exp(-u * Lg)) - 1) \
        * sp.exp(u * Lg) / g
    cc = mp.cos(alpha)
    out = []
    for k in range(K + 1):
        f = sp.lambdify((g, Lg, Lp, c), expr.subs(u, 0), "mpmath")
        val = mp.quad(lambda x: f(x, mp.log(x), mp.log1p(x), cc), [0, mp.mpf(1) / 64, mp.mpf(1) / 8, mp.mpf(1) / 2, 1])
        out.append(val / mp.factorial(k))
        expr = sp.diff(expr, u)
    return out


def sec_corner():
    print("\n== 8. Corner coefficients at interior angle alpha (dps 70) ==")
    mp.mp.dps = 70
    C = consts()
    pi = mp.pi
    l2, l3, pi2, z3 = mp.log(2), mp.log(3), mp.pi**2, mp.zeta(3)
    Li2_13 = mp.polylog(2, mp.mpf(1) / 3)
    # 50-digit values from an earlier independent computation (interior angles)
    PREV = {"3pi/4": "0.07350854650766740327318436636766157962413565493423",
          "2pi/3": "0.12579567893797323375386819768385949481789421068628",
          "pi/2": "0.25261909620704479461836352096565486053362715072438",
          "pi/3": "0.37802988488257911145014904967078533671288591101805",
          "pi/4": "0.42748316191642018108039374021131783899370796572390",
          "pi/6": "0.46196590551845269597956607891396239753617124296662"}
    angles = {"3pi/4": 3 * pi / 4, "2pi/3": 2 * pi / 3, "pi/2": pi / 2, "pi/3": pi / 3, "pi/4": pi / 4, "pi/6": pi / 6}
    res = {}
    for name, al in angles.items():
        t = time.time()
        Mk = corner_coeffs(al, 2)
        res[name] = Mk
        cA = mp.cos(al / 2)
        Kal = -(1 + cA) * mp.log(1 + cA) - (1 - cA) * mp.log(1 - cA)
        print(f"  alpha = {name}: M0 = {mp.nstr(Mk[0], 25)}, M1 = {mp.nstr(Mk[1], 25)}, M2 = {mp.nstr(Mk[2], 25)}"
              f"  [{time.time() - t:.0f}s]")
        check(f"M0(alpha) M0({name}) = -(1+c)log(1+c)-(1-c)log(1-c)", Mk[0], Kal, 60)
        check(f"Corner M1({name}) vs earlier 50-digit value", Mk[1], mp.mpf(PREV[name]), 48)
    check("Corner M1(pi/2) = Theorem 1", res["pi/2"][1], M1_closed(C), 60)
    check("Corner M2(pi/2) = Theorem 2", res["pi/2"][2], M2_closed(C), 60)
    rhs_a = (pi2 + 48 * l2**2 - 6 * l3**2 - 24 * l2 * l3 - 12 * Li2_13) / 24
    check("Corner 24 M1(2pi/3) = pi^2+48l2^2-6l3^2-24l2l3-12Li2(1/3)", res["2pi/3"][1], rhs_a, 60)
    s3 = mp.sqrt(3)
    L3 = mp.log(2 + s3)
    rhs_d = (8 * l2**2 - 2 * l3**2 - 4 * L3**2 - 4 * Li2_13
             - s3 * (2 * pi2 + l2**2 + l3**2 - 5 * L3**2 - 2 * l2 * l3 + 6 * l2 * L3 - 2 * l3 * L3 + 2 * Li2_13
                     - 24 * mp.polylog(2, s3 - 1) + 8 * mp.polylog(2, (s3 - 1) / 2))) / 4
    check("Corner 4 M1(pi/3) = [Q(sqrt3) form]", res["pi/3"][1], rhs_d, 60)
    Li3_13, Li3_m12 = mp.polylog(3, mp.mpf(1) / 3), mp.polylog(3, -mp.mpf(1) / 2)
    br = 48 * (l2**2 / 2 - mp.polylog(2, -mp.mpf(1) / 2))
    rhs_b = (9 * z3 - 20 * pi2 * l2 + 12 * pi2 * l3 - 4 * l2**3 + 60 * l2**2 * l3 - 84 * l2 * l3**2 + 30 * l3**3
             - 168 * l2 * Li2_13 + 180 * l3 * Li2_13 + 360 * Li3_13 + 408 * Li3_m12 + br) / 48
    check("Corner 48 M2(2pi/3) = [Q form with zeta3, Li3(1/3), Li3(-1/2)]", res["2pi/3"][2], rhs_b, 60)
    check("Corner bracket identity 48l2^2+24l3^2-48l2l3+48Li2(1/3) = 48(l2^2/2-Li2(-1/2))",
          48 * l2**2 + 24 * l3**2 - 48 * l2 * l3 + 48 * Li2_13, br, 60)

    # ---- quartic cases 3pi/4 and pi/6 (closed forms transcribed from Table II of the paper)
    Q = mp.mpf
    Li2 = lambda x: mp.polylog(2, x)
    z2 = mp.zeta(2)
    s2, s3, s6 = mp.sqrt(2), mp.sqrt(3), mp.sqrt(6)
    th, ph = mp.sqrt(2 - s2), mp.sqrt(2 + s2)
    li_34 = [(-1 + th / 2, 5 - 6 * th - 2 * th**2 + 2 * th**3), (-1 + th / 2, 1 - th - th**2 / 2 + th**3 / 2),
             (1 + th / 2, -3 + 4 * th + th**2 - th**3), (1 + th / 2, Q(1) / 2 - th / 2 + th**3 / 4),
             (2 - th, 1 - 2 * th**2 + th**3), (-2 * th, 2 - 2 * th - th**2 + th**3), (-2 - th, th / 2),
             (-th, 2 - 2 * th - th**2 / 2 + th**3 / 2), (-2 * th, -1 + 4 * th**2 - 2 * th**3),
             (2, -7 + 10 * th + 2 * th**2 - 3 * th**3)]
    l = [mp.log(th), mp.log(1 + s2), mp.log(1 + th), mp.log(1 + ph)]
    pr_34 = {(0, 0): Q(39) / 2 + 3 * th / 4, (0, 1): 11 + 31 * th / 2, (0, 2): -6 + 9 * th, (0, 3): 3 - 29 * th / 2,
             (1, 1): 2 + th, (1, 2): -5 - 9 * th / 2, (1, 3): 2 * th, (2, 2): Q(1) / 2 - 11 * th / 4,
             (2, 3): 3 + 11 * th / 2, (3, 3): Q(1) / 2 - 11 * th / 4}
    M1_34 = (mp.fsum(cf * Li2(x) for cf, x in li_34) + (5 * th / 2 - 1) * z2
             + mp.fsum(cf * l[i] * l[j] for (i, j), cf in pr_34.items()))
    check("TableII M1(3pi/4) = closed form (Q(cos 3pi/8))", res["3pi/4"][1], M1_34, 65)
    c34 = mp.cos(3 * pi / 8)
    check("TableII   letter 1-2th^2+th^3 = c/(1+c), c = cos(3pi/8)", 1 - 2 * th**2 + th**3, c34 / (1 + c34), 65)
    check("TableII   letter -1+4th^2-2th^3 = tan^2(3pi/16)", -1 + 4 * th**2 - 2 * th**3, mp.tan(3 * pi / 16)**2, 65)
    check("TableII   letter -7+10th+2th^2-3th^3 = phi/(2+phi)", -7 + 10 * th + 2 * th**2 - 3 * th**3, ph / (2 + ph), 65)
    sg = s2 + s6
    li_6 = [(-sg, 15 - 10 * s2 + 8 * s3 - 6 * s6), (2 + sg / 2, 1 - sg / 4), (-1 + sg / 4, 5 - 2 * s6),
            (-sg / 2, Q(3) / 2 + s2 / 2 - s3 / 2 - s6 / 2), (sg, -2 + s2 - s3 + s6), (-1 - sg / 4, Q(1) / 2 - s2 / 4),
            (2, -7 - 5 * s2 + 4 * s3 + 3 * s6), (-1 - sg / 4, 3 + s2 - s3 - s6),
            (-1 + sg / 4, Q(3) / 2 - s2 / 2 + s3 / 2 - s6 / 2), (2 - sg / 2, -7 + 5 * s2 - 4 * s3 + 3 * s6)]
    al6 = mp.sqrt(2 + s3) - 1
    m = [mp.log(al6), mp.log(1 + s2), mp.log(2 + s3), mp.log(s2 + s3)]
    pr_6 = {(0, 0): 28 - sg, (0, 1): 2 + 13 * sg / 2, (0, 2): -9 - 5 * sg / 4, (0, 3): 30 + 7 * sg / 2,
            (1, 1): Q(3) / 2 - 9 * sg / 8, (1, 2): -Q(7) / 2 - 3 * sg / 8, (1, 3): -1 + 7 * sg / 4,
            (2, 2): Q(9) / 8 + 5 * sg / 32, (2, 3): -Q(11) / 2 - 7 * sg / 8, (3, 3): Q(15) / 2 + 11 * sg / 8}
    c6 = mp.cos(pi / 12)
    M1_6 = (mp.fsum(cf * Li2(x) for cf, x in li_6) + (c6 - 1) * z2
            + mp.fsum(cf * m[i] * m[j] for (i, j), cf in pr_6.items()))
    check("TableII M1(pi/6) = closed form (Q(sqrt2,sqrt3))", res["pi/6"][1], M1_6, 65)
    check("TableII   c = cos(pi/12) = (sqrt6+sqrt2)/4; alpha = 2c-1", al6, 2 * c6 - 1, 65)
    check("TableII   letter -2+s2-s3+s6 = tan(pi/24)", -2 + s2 - s3 + s6, mp.tan(pi / 24), 65)
    check("TableII   letter -7+5s2-4s3+3s6 = c/(1+c)", -7 + 5 * s2 - 4 * s3 + 3 * s6, c6 / (1 + c6), 65)
    t1 = 1 + 2 * mp.sin(pi / 12)
    check("TableII   letter -7-5s2+4s3+3s6 = (t1-1)/(t1+1), t1 = 1+2sin(pi/12)", -7 - 5 * s2 + 4 * s3 + 3 * s6,
          (t1 - 1) / (t1 + 1), 65)
    # M2, M3 at the quartic angles vs earlier 61-digit values
    d10b = {"3pi/4": ("-0.0053731015427211737654998917149816026461899163711776354058878",
                      "0.0060901742186917557579690789800261569379082184706542129260597"),
            "pi/6": ("-0.0761405079464313619726596784453231945776006850073666905590217",
                     "0.0351375430991826292429345021274375869660672136078787624622955")}
    for name, (v2, v3) in d10b.items():
        Mk = corner_coeffs(angles[name], 3)
        check(f"TableII M2({name}) vs printed value", Mk[2], mp.mpf(v2), 58)
        check(f"TableII M3({name}) vs printed value", Mk[3], mp.mpf(v3), 58)


# -------------------------------------------------------------------------------------------
# Section 9: 1200-digit checks of Theorems 1, 2 and Conjecture 4 against the note's stored values
# -------------------------------------------------------------------------------------------
NOTE_OUT = HERE / "squircle_square_end" / "outputs"


def sec_hi():
    print("\n== 9. Theorems 1-4 vs the stored 1200-digit M_k (methods A, B; dps 1215) ==")
    mp.mp.dps = 1215
    t = time.time()
    C = consts()
    st = {}
    for k in range(4):
        for m in "AB":
            s = open(f"{NOTE_OUT}/M{k}_method{m}_1200.txt").read().strip()     # read only
            st[(k, m)] = mp.mpf(s)
    closed = {0: M0_closed(C), 1: M1_closed(C), 2: M2_closed(C)}
    for k in range(3):
        for m in "AB":
            check(f"F{2 + k}  M{k} closed form vs stored 1200-digit method {m}", st[(k, m)], closed[k], 1190)
    check("Thm3  M3: stored method A vs method B", st[(3, "A")], st[(3, "B")], 1190)
    R, R0, R1 = conj4(C)
    W3p, T = W3p_closed(C), T_closed(C)
    for m in "AB":
        W4 = st[(3, m)] - W3p - T
        check(f"Thm4  Conjecture 4 RHS vs W4 = M3(method {m}) - W3' - T", W4, R, 1190)
    show("(W4 - R)/pi^4 with M3 from method A", (st[(3, "A")] - W3p - T - R) / mp.pi**4, 5)
    print(f"  [{time.time() - t:.0f}s]")


# -------------------------------------------------------------------------------------------
# Section 10: re-run (NOT independent) of the certified ball-arithmetic bound for c0, c1
# -------------------------------------------------------------------------------------------
def sec_cert():
    import re
    import subprocess
    print("\n== 10. Re-run of weight4_M3/w4_certified.py 1100 1200 (python-flint balls; RE-RUN, not independent) ==")
    script = str(HERE / "weight4_M3" / "w4_certified.py")
    t = time.time()
    out = subprocess.run([sys.executable, "-B", script, "1100", "1200"], capture_output=True, text=True,
                         cwd=str(HERE / "weight4_M3"), timeout=1800)
    txt = out.stdout
    for line in txt.splitlines():
        if line.startswith(("PREC", "c0 =", "c1 =")) or "CERTIFIED" in line:
            print("   | " + line.strip())
    print(f"  [returncode {out.returncode}, {time.time() - t:.0f}s]")
    for name, k in (("c0", 322), ("c1", 321)):
        ok = f"CERTIFIED |{name}| < 10^-{k} : True" in txt
        LOG.append((f"Thm4  certified |{name}| < 1e-{k} (re-run)", k if ok else -1, k, ok))
        print(f"  F6  certified |{name}| < 1e-{k} (re-run)  {'PASS' if ok else '**FAIL**'}")
    rad = {}
    for name in ("c0", "c1"):
        mm = re.search(name + r" = \[\+/- ([0-9.e+-]+)\]", txt)
        rad[name] = mp.mpf(mm.group(1)) if mm else mp.inf
    mp.mp.dps = 30
    bound = mp.pi**4 * (rad["c0"] + mp.sqrt(2) * rad["c1"])
    print(f"  => |W4 - R| = pi^4 |c0 + sqrt2 c1| <= {mp.nstr(bound, 4)}  (paper: < 1e-319)")
    LOG.append(("Thm4  |W4 - R| < 1e-319 from the printed balls", 0, 0, bound < mp.mpf(10)**-319))


# -------------------------------------------------------------------------------------------
# Section 11: pi_p (l^p length of the l^p unit circle, halved): expansion at p = infinity 
# -------------------------------------------------------------------------------------------
def _smul(a, b, K):
    return [mp.fsum(a[i] * b[n - i] for i in range(n + 1)) for n in range(K + 1)]


def _sexp(a, K):
    f = [mp.exp(a[0])] + [mp.mpf(0)] * K
    for n in range(1, K + 1):
        f[n] = mp.fsum(k * a[k] * f[n - k] for k in range(1, n + 1)) / n
    return f


def _slog(a, K):
    g = [mp.log(a[0])] + [mp.mpf(0)] * K
    for n in range(1, K + 1):
        g[n] = (n * a[n] - mp.fsum(k * g[k] * a[n - k] for k in range(1, n))) / (n * a[0])
    return g


def pip_exact(p):
    """pi_p from the curve itself: pi_p = 4 int_0^{2^(-1/p)} (1 + x^(p(p-1)) (1-x^p)^(1-p))^(1/p) dx."""
    p = mp.mpf(p)
    f = lambda x: (1 + x**(p * (p - 1)) * (1 - x**p)**(1 - p))**(1 / p)
    return 4 * mp.quad(f, [0, 2**(-1 / p) / 2, 2**(-1 / p)])


def pip_R(eps):
    """Own form (derived here): pi_p = 4*2^(-eps) + 4 eps^2 int_0^inf G(eps,u) du, eps = 1/p,
    G = e^{eps u} (1+e^{eps u})^(-1-eps) [ (1+e^{-u(1-eps)})^eps - 1 ]."""
    e = mp.mpf(eps)
    G = lambda u: mp.exp(e * u) * (1 + mp.exp(e * u))**(-1 - e) * (mp.expm1(e * mp.log1p(mp.exp(-u * (1 - e)))))
    return 4 * 2**(-e) + 4 * e**2 * mp.quad(G, [0, 1, 4, 16, 64, 256, mp.inf])


def sec_pip(K=26):
    print(f"\n== 11. pi_p: expansion coefficients N_0..N_{K} by series arithmetic in eps (own form, dps 60) ==")
    mp.mp.dps = 60
    lam, pi2, z3, z5 = mp.log(2), mp.pi**2, mp.zeta(3), mp.zeta(5)
    check("Thm8 own form = curve definition, p = 2 (pi_2 = pi)", pip_R(mp.mpf(1) / 2), mp.pi, 55)
    check("Thm8 own form = curve definition, p = 7", pip_R(mp.mpf(1) / 7), pip_exact(7), 50)
    check("Thm8 duality pi_3 = pi_{3/2} (curve definition)", pip_exact(3), pip_exact(mp.mpf(3) / 2), 50)
    check("Thm8 duality pi_5 = pi_{5/4} (curve definition)", pip_exact(5), pip_exact(mp.mpf(5) / 4), 50)
    # Taylor coefficients in eps of G(eps,u), per node u (cached), then integrate each coefficient
    E = [1 / mp.factorial(j) for j in range(K + 1)]
    lser = _slog([mp.mpf(2)] + E[1:], K)            # log(1+e^x) = sum lser_j x^j
    cache = {}

    def gser(u):
        if u in cache:
            return cache[u]
        Lu = [lser[j] * u**j for j in range(K + 1)]
        A = [-Lu[0]] + [(u if j == 1 else 0) - Lu[j] - Lu[j - 1] for j in range(1, K + 1)]
        w = mp.exp(-u)
        B = _slog([1 + w] + [w * u**j * E[j] for j in range(1, K + 1)], K)
        br = _sexp([mp.mpf(0)] + B[:K], K)
        br[0] -= 1
        g = _smul(_sexp(A, K), br, K)
        cache[u] = g
        return g
    t = time.time()
    N = [-lam]
    for k in range(1, K + 1):
        Ik = mp.quad(lambda u: gser(u)[k - 1], [0, 1, 4, 16, 64, 256, mp.inf])
        N.append((-lam)**(k + 1) / mp.factorial(k + 1) + Ik)
    print(f"  [N_0..N_{K} in {time.time() - t:.0f}s, {len(cache)} nodes]")
    closed = {0: -lam, 1: lam**2 / 2, 2: pi2 / 24 - lam**3 / 6,
              3: lam**4 / 24 + pi2 / 24 * (1 - lam) + z3 / 4,
              4: -lam**5 / 120 + pi2 / 24 * (1 - lam + lam**2 / 2) + z3 / 4 * (1 - lam) + pi2**2 / 2880,
              5: (lam**6 / 720 + pi2 / 24 * (1 - lam + lam**2 / 2 - lam**3 / 6) + z3 / 4 * (1 - lam + lam**2 / 2)
                  + pi2**2 / 2880 * (1 - lam) - 7 * pi2**2 / 1920 - z3 * pi2 / 48 + z5 / 8)}
    for k, v in closed.items():
        check(f"Thm8 N_{k} closed form  vs own series coefficient", N[k], v, 50)
    d12 = {6: "-0.9326682134840535164505582815863254270674", 7: "-1.472571756898086547743541305319946712566",
           8: "0.1221389382135391530907648898935863780245", 9: "4.700293562406214840901019828506766032978",
           10: "4.4249", 11: "-18.181", 12: "-30.262", 13: "149.34", 14: "299.17", 15: "-2073.5",
           20: "-2.5866e6", 25: "7.3762e10"}
    for k, v in d12.items():
        nd = len(v.split("e")[0].replace("-", "").replace(".", "").lstrip("0"))
        check(f"Thm8 N_{k} vs stored table", N[k], mp.mpf(v), nd - 2)
    print("   k    N_k pi^(k-1) / (2 (k-2)!)   (stored: 1.72897, 1.40548, 1.28331, 1.21873 at k = 13, 17, 21, 25)")
    for k in [9, 13, 17, 21, 25]:
        print(f"  {k:3d}  {mp.nstr(N[k] * mp.pi**(k - 1) / (2 * mp.factorial(k - 2)), 8)}")
    check("Thm8 growth: N_25 pi^24/(2*23!) vs stored 1.21873", N[25] * mp.pi**24 / (2 * mp.factorial(23)),
          mp.mpf("1.21873"), 5)
    # asymptotic check against the exact value at large p:  (pi_p/4 - sum_{k<=8} N_k p^-(k+1)) p^10 -> N_9
    for p in [200, 400]:
        e = mp.mpf(1) / p
        rem = (pip_R(e) / 4 - 1 - mp.fsum(N[k] * e**(k + 1) for k in range(9))) / e**10
        print(f"  p = {p}: (pi_p/4 - series through N_8) * p^10 = {mp.nstr(rem, 10)}   (N_9 = {mp.nstr(N[9], 10)})")
    e = mp.mpf(1) / 400
    rem = (pip_R(e) / 4 - 1 - mp.fsum(N[k] * e**(k + 1) for k in range(9))) / e**10
    check("Thm8 exact pi_400 vs series through N_8: remainder*p^10 ~ N_9", rem, N[9], 1)


# -------------------------------------------------------------------------------------------
# Section 12: superellipsoid surface representation (D4) at p = 1, 2
# -------------------------------------------------------------------------------------------
def sec_ellipsoid():
    print("\n== 12. Superellipsoid |x|^p+|y|^p+|z|^p = 1: representation S/8 = int_simplex ... (dps 20) ==")
    mp.mp.dps = 20

    def S(p):
        e = 1 / mp.mpf(p)

        def inner(A):
            w = 1 - A
            def f(v):                      # B = w v, T = w (1 - v): no rounding at the hypotenuse
                B, T = w * v, w * (1 - v)
                return w * e * e * (A * B * T)**(e - 1) * mp.sqrt(A**(2 - 2 * e) + B**(2 - 2 * e) + T**(2 - 2 * e))
            return mp.quad(f, [0, mp.mpf(1) / 2, 1])
        return 8 * mp.quad(inner, [0, mp.mpf(1) / 2, 1])
    check("Thm14 representation at p = 1: S = 4 sqrt3 (octahedron)", S(1), 4 * mp.sqrt(3), 15)
    check("Thm14 representation at p = 2: S = 4 pi (sphere)", S(2), 4 * mp.pi, 8)


# -------------------------------------------------------------------------------------------
# Section 13: M_1(alpha) at every interior angle (Theorem 11) vs direct quadrature
# -------------------------------------------------------------------------------------------
def sec_cornerthm():
    print("\n== 13. M_1(alpha) for all angles: real-dilogarithm theorem vs quadrature (dps 40) ==")
    mp.mp.dps = 40
    L2 = lambda x: mp.re(mp.polylog(2, x))
    lg = mp.log

    def M1_thm(a):
        c, s = mp.cos(a / 2), mp.sin(a / 2)
        Xm, Xp = (1 + s - c) / (1 - c), (1 + s + c) / (1 + c)
        U = (-L2(-1 / s) + L2(-s / (1 - c)) + L2(-s / (1 + c)) - L2(-(1 + s) / c) - L2(-s * (1 + s) / c**2)
             + L2(-s * (1 + s) / (c * (1 + c))) + c * L2(-(1 + s - c) / c) + c * L2(-(1 + s - c) / (c * (1 - c)))
             - 2 * c * L2(-(1 + s - c) / (2 * c)) - 2 * L2(1 / (1 + s)) + L2(c**2 / (1 + s)) + L2(c / (1 + s))
             + L2(s / (1 + s)) - L2(c * (1 - c) / (s * (1 + s))) - c * L2((1 - c) / (1 + s - c))
             + c * L2((1 + c) / (1 + s + c)) + c * L2(c / (1 + s + c)) + c * L2(c * (1 + c) / (1 + s + c))
             - 2 * c * L2(2 * c / (1 + s + c)))
        V = ((1 - c) * L2(-1 / c) + L2(1 - 1 / c**2) - (1 + c) * L2(-(1 - c) / c) + 2 * c * L2(-(1 - c) / (2 * c))
             - L2(c**2) - (1 + c) * L2(c) + (1 - c) * L2(c / (1 + c)) + 2 * c * L2(2 * c / (1 + c)))
        Lp = (lg(1 + s) * lg(c / (1 + s)) + lg(c / (1 + c)) * lg((1 + s) / s) + c * lg(4) * lg(Xm)
              + c / 2 * lg(Xm)**2 - c / 2 * lg(Xp) * lg((1 + c) * (1 + s + c)))
        return mp.zeta(2) + U + V + Lp

    def M1_quad(a):                       # direct quadrature of the u-derivative at u = 0
        ca = mp.cos(a)

        def f(g):
            S0 = 1 - 2 * ca * g + g * g
            f0 = mp.sqrt(S0) / ((1 + g) * g)
            l = mp.log(g)
            return f0 * (-mp.log(1 + g) + l + (-l * g * g + ca * l * g) / S0) - l / g
        return mp.quad(f, [0, mp.mpf(1) / 4, mp.mpf(1) / 2, 1])

    for a in [mp.mpf("0.05"), mp.pi / 5, mp.pi / 2, 2 * mp.pi / 3, mp.mpf("2.5"), mp.mpf("3.1")]:
        check(f"Thm11 M1(alpha) theorem vs quadrature, alpha = {mp.nstr(a, 6)}", M1_thm(a), M1_quad(a), 35)


# -------------------------------------------------------------------------------------------
# Section 14: smoothed regular 2m-gons at p = 1 (Theorem 12, Observation 2)
# -------------------------------------------------------------------------------------------
def sec_hexagon():
    print("\n== 14. Smoothed 2m-gons at p = 1: P(1), the P'(1) law, hexagon derivatives (dps 50) ==")
    mp.mp.dps = 50

    def perim(m, p):                      # polar perimeter of {sum_i |<n_i,x>|^p <= 1}, n_i = e(i pi/m)
        ang = [i * mp.pi / m for i in range(1, m + 1)]

        def f(th):
            c = [mp.cos(th - a) for a in ang]
            s = [mp.sin(th - a) for a in ang]
            S = mp.fsum(abs(x)**p for x in c)
            Sp = mp.fsum(-p * abs(x)**(p - 1) * mp.sign(x) * y for x, y in zip(c, s))
            rho = S**(-1 / p)
            drho = -(1 / p) * S**(-1 / p - 1) * Sp
            return mp.sqrt(rho**2 + drho**2)
        kinks = sorted(set([mp.mpf(0), mp.pi] + [(a + mp.pi / 2) % mp.pi for a in ang]))
        return 2 * mp.quad(f, kinks)

    for m in (2, 3, 4):
        check(f"Thm12 P(1) = 4m(sec(pi/2m) - cos(pi/2m)), m = {m}", perim(m, 1),
              4 * m * (mp.sec(mp.pi / (2 * m)) - mp.cos(mp.pi / (2 * m))), 40)
    m, h, K = 3, mp.mpf("0.004"), 7
    xs = [k * h for k in range(-K, K + 1)]
    ys = [perim(m, 1 + x) for x in xs]
    co = mp.lu_solve(mp.matrix([[x**j for j in range(2 * K + 1)] for x in xs]), mp.matrix(ys))
    d1, d2, d3 = co[1], 2 * co[2], 6 * co[3]
    r3, L = mp.sqrt(3), mp.log(2)
    V = mp.cot(mp.pi / 6)
    R1 = mp.fsum(mp.sin(i * mp.pi / 3) * mp.log(V / mp.sin(i * mp.pi / 3)) for i in (1, 2)) / V
    check("Thm12 hexagon: R_1 = log 2", R1, L, 45)
    check("Thm12 hexagon: P'(1) = P(1) R_1 = 2 sqrt3 log2 (finite differences)", d1, 2 * r3 * R1, 17)
    check("Obs2 hexagon P''(1) closed form (finite differences)", d2,
          r3 / 8 * mp.pi**2 + 2 * r3 * L**2 - 4 * r3 * L, 15)
    check("Obs2 hexagon P'''(1) closed form (finite differences)", d3,
          r3 / 8 * (16 * L**3 + 3 * mp.pi**2 * L - 96 * L**2 + 96 * L - 9 * mp.zeta(3) - 6 * mp.pi**2), 12)

SECTIONS = {"squircle": sec_squircle, "vform": sec_vform, "lemmas": sec_lemmas, "rhombus": sec_rhombus,
            "borel": sec_borel, "mk": sec_mk, "area": sec_area, "corner": sec_corner, "hi": sec_hi,
            "cert": sec_cert, "pip": sec_pip, "ellipsoid": sec_ellipsoid,
            "cornerthm": sec_cornerthm, "hexagon": sec_hexagon}

if __name__ == "__main__":
    which = sys.argv[1:] or list(SECTIONS)
    M = None
    for s in which:
        if s == "squircle":
            M = sec_squircle()
        elif s == "vform":
            sec_vform(M)
        else:
            SECTIONS[s]()
    nfail = sum(1 for _, _, _, ok in LOG if not ok)
    print(f"\nSUMMARY: {len(LOG)} checks, {nfail} failed, {time.time() - T0:.0f}s total")
    for tag, dg, need, ok in LOG:
        if not ok:
            print(f"  FAILED: {tag}  (digits {dg}, need {need})")
