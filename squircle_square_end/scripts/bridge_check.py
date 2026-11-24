"""Remark (L(chi_8,3)) and Appendix B of the note.

(1) Cl3(pi/4) = (sqrt2/2) L(chi8,3) - 3 zeta(3)/256   (elementary: split n mod 8)  -- numerical check
(2) 48 chi3(sqrt2-1) = 32 sqrt2 L(chi8,3) - 3 pi^2 L + 4 L^3                        -- numerical check ONLY
(3) symbolic: given (1),(2), Theorem 2 (real-Li3 form) == L(chi8,3) form == Cl3(pi/4) form.
(4) the proved closed form vs the 1200-digit value of M2 (if available).
"""
import os, sys
import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "outputs")
DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 1200
mp.mp.dps = DPS + 20
lines = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


r2 = mp.sqrt(2); b = r2 - 1; L = mp.log(1 + r2); lam = mp.log(2); pi2 = mp.pi ** 2; z3 = mp.zeta(3)
L8 = mp.dirichlet(3, [0, 1, 0, -1, 0, -1, 0, 1])
Cl3 = mp.clcos(3, mp.pi / 4)
Li3b, Li3mb, Li3h, Li2b = mp.polylog(3, b), mp.polylog(3, -b), mp.polylog(3, 1 / r2), mp.polylog(2, b)
chi3 = (Li3b - Li3mb) / 2
say(f"working precision {DPS} digits")
say("(1) Cl3(pi/4) - (sqrt2/2 L8 - 3 z3/256)       =", mp.nstr(Cl3 - (r2 / 2 * L8 - 3 * z3 / 256), 5))
say("(2) 48 chi3(b) - (32 sqrt2 L8 - 3 pi^2 L + 4L^3) =", mp.nstr(48 * chi3 - (32 * r2 * L8 - 3 * pi2 * L + 4 * L ** 3), 5))

S2 = sp.sqrt(2)
Z, LI3B, LI3MB, LI3H, LI2B, LL, LAM, PI, L8s, CL3 = sp.symbols("z3 Li3b Li3mb Li3h Li2b L lam pi L8 Cl3")
R = sp.Rational
thm2 = ((PI ** 2 + 6 * LAM ** 2) / 24 + (R(189, 64) * S2 - R(19, 16)) * Z + (11 + 15 * S2) * LI3B + (11 - 15 * S2) * LI3MB
        - 27 * S2 * LI3H + (12 * LL - R(7, 2) * S2 * LAM) * LI2B + (R(15, 8) * S2 - R(7, 12)) * PI ** 2 * LL
        - (R(1, 24) + R(65, 48) * S2) * PI ** 2 * LAM + (R(1, 6) + R(9, 16) * S2) * LAM ** 3 + R(11, 8) * S2 * LAM ** 2 * LL
        - (1 + R(7, 4) * S2) * LAM * LL ** 2 + (R(25, 6) - R(7, 3) * S2) * LL ** 3)
formB = (PI ** 2 + 6 * LAM ** 2) / 24 + (-228 * Z + 4224 * LI3B + 7680 * L8s - 8 * PI ** 2 * LAM + 152 * PI ** 2 * LL
         + 32 * LAM ** 3 - 192 * LAM * LL ** 2 + 448 * LL ** 3 + 2304 * LI2B * LL
         + S2 * (567 * Z - 5184 * LI3H - 2816 * L8s - 260 * PI ** 2 * LAM + 108 * LAM ** 3 + 264 * LAM ** 2 * LL
                 - 336 * LAM * LL ** 2 + 32 * LL ** 3 - 672 * LI2B * LAM)) / 192
formA = (PI ** 2 / 24 + LAM ** 2 / 4 + LAM ** 3 / 6 - LAM * LL ** 2 + R(7, 3) * LL ** 3 - PI ** 2 * LAM / 24
         + R(19, 24) * PI ** 2 * LL + 12 * LL * LI2B - R(49, 32) * Z - R(88, 3) * CL3 + 22 * LI3B
         + S2 * (R(9, 16) * LAM ** 3 + R(11, 8) * LAM ** 2 * LL - R(7, 4) * LAM * LL ** 2 + LL ** 3 / 6
                 - R(65, 48) * PI ** 2 * LAM - R(7, 2) * LAM * LI2B + R(219, 64) * Z + 40 * CL3 - 27 * LI3H))
chi3_sym = (32 * S2 * L8s - 3 * PI ** 2 * LL + 4 * LL ** 3) / 48
thm2_L8 = sp.expand(thm2.subs(LI3MB, LI3B - 2 * chi3_sym))
say("(3a) Theorem 2 [Li3(-b) -> Li3(b) - 2 chi3, (2)]  -  L(chi8,3) form =", sp.simplify(sp.expand(thm2_L8 - formB)))
formA_L8 = sp.expand(formA.subs(CL3, S2 / 2 * L8s - R(3, 256) * Z))
say("(3b) Cl3 form [Cl3 -> (1)]  -  L(chi8,3) form                        =", sp.simplify(sp.expand(formA_L8 - formB)))
say("     coefficient of L(chi8,3) in the L-form:", sp.radsimp(sp.expand(formB).coeff(L8s)))

num = sp.lambdify((Z, LI3B, LI3MB, LI3H, LI2B, LL, LAM, PI), thm2, modules="mpmath")
M2_cf = num(z3, Li3b, Li3mb, Li3h, Li2b, L, lam, mp.pi)
for meth in ("A", "B"):
    f = os.path.join(OUT, f"M2_method{meth}_1200.txt")
    if os.path.exists(f):
        ref = mp.mpf(open(f).read().strip())
        d = M2_cf - ref
        say(f"(4) Theorem 2 closed form - M2 (method {meth}, 1200 digits) = {mp.nstr(d, 5)}"
            f"  -> {int(-mp.log10(abs(d / ref))) if d else 'all'} digits")
with mp.workdps(60):
    say("    M2 =", mp.nstr(M2_cf, 55))
with open(os.path.join(OUT, "bridge_check.out"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
