"""Assemble the experimental closed form of M3 = T + W3' + W4 (T, W3' proved; W4 from the LLL relation of
m3_closed_form.py), check it against the 1200-digit M3 from methods A and B, and emit LaTeX.
"""
import os, sys, pickle
import mpmath as mp
import polylog_exact as PE   # noqa: F401 (import before setting precision)
import sympy as sp
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
OUTD = os.path.join(HERE, "..", "outputs")
mp.mp.dps = 1220
rel = pickle.load(open(os.path.join(OUTD, "m3_relation.pkl"), "rb"))
tails = pickle.load(open(os.path.join(OUTD, "m3_tails_exact.pkl"), "rb"))
mp.mp.dps = 1220
S2 = sp.sqrt(2)
sy = {k: sp.Symbol(k) for k in ["lam", "L", "pi2", "Li2b", "z3", "Li3b", "Li3mb", "Li3h"]}


def qs(c):
    return sp.Rational(c.p.numerator, c.p.denominator) + sp.Rational(c.q.numerator, c.q.denominator) * S2


def P_to_sympy(p):
    from qsqrt2 import B, QS2
    out = 0
    for k, c in p.d.items():
        term = qs(c)
        for a in k:
            if isinstance(a, str):
                term *= sy[a]
            elif a[0] == 'Li2':
                assert a[1] == B; term *= sy["Li2b"]
            else:
                z = a[1]
                key = {B: "Li3b", -B: "Li3mb", QS2(0, Fr(1, 2)): "Li3h"}[z]
                term *= sy[key]
        out += term
    return sp.expand(out)


T = P_to_sympy(tails["T"]); W3p = P_to_sympy(tails["W3p"])
vec, kept, prods = rel["vec"], rel["kept"], rel["prods"]
c0 = vec[0]
W4 = 0
idx = 1
for p in prods:
    mono = sp.Mul(*[sy[x] for x in p])
    W4 += sp.Rational(-vec[idx], c0) * mono + sp.Rational(-vec[idx + 1], c0) * S2 * mono
    idx += 2
li4s = {}
for z in kept:
    s_ = sp.Symbol(f"Li4_{len(li4s)}")
    li4s[z] = s_
    W4 += sp.Rational(-vec[idx], c0) * s_ + sp.Rational(-vec[idx + 1], c0) * S2 * s_
    idx += 2
M3 = sp.expand(T + W3p + W4)

r2 = mp.sqrt(2); b = r2 - 1
vals = {"lam": mp.log(2), "L": mp.log(1 + r2), "pi2": mp.pi ** 2, "Li2b": mp.polylog(2, b), "z3": mp.zeta(3),
        "Li3b": mp.polylog(3, b), "Li3mb": mp.polylog(3, -b), "Li3h": mp.polylog(3, 1 / r2)}
syms = list(sy.values()) + list(li4s.values())
f = sp.lambdify(syms, M3, modules="mpmath")
num = f(*[vals[k] for k in sy], *[mp.polylog(4, PE.DE.float_exact(z)) for z in kept])
lines = []
for m in "AB":
    ref = mp.mpf(open(os.path.join(OUTD, f"M3_method{m}_1200.txt")).read().strip())
    d = num - ref
    lines.append(f"closed form - M3(method {m}, 1200 digits) = {mp.nstr(d, 3)}")
with mp.workdps(70):
    lines.append("M3 (closed form) = " + mp.nstr(num, 65))

# LaTeX: group by constant
names = {"lam": r"\lambda", "L": "L", "pi2": r"\pi^2", "Li2b": r"\Li_2(b)", "z3": r"\zeta(3)",
         "Li3b": r"\Li_3(b)", "Li3mb": r"\Li_3(-b)", "Li3h": r"\Li_3(\tfrac1{\sqrt2})"}
from qsqrt2 import QS2, B, SQ2
zname = {QS2(Fr(1, 2)): r"\tfrac12", B: "b", -B: "-b", QS2(0, Fr(1, 2)): r"\tfrac1{\sqrt2}", -B * B: "-b^2",
         QS2(2) - SQ2: r"2-\sqrt2", QS2(1) - QS2(0, Fr(1, 2)): r"1-\tfrac1{\sqrt2}",
         QS2(-1) + QS2(0, Fr(1, 2)): r"\tfrac1{\sqrt2}-1"}
for z, s_ in li4s.items():
    names[str(s_)] = r"\Li_4(" + zname[z] + ")"


def coef_tex(c):
    c = sp.nsimplify(c, [S2]) if False else sp.radsimp(c)
    p = sp.Rational(c.subs(S2, 0)); q = sp.Rational(sp.expand((c - p) / S2))
    def fr(x):
        x = sp.Rational(x)
        return (str(abs(x.p)) if x.q == 1 else r"\tfrac{%d}{%d}" % (abs(x.p), x.q))
    parts = []
    if p != 0:
        parts.append(("-" if p < 0 else "+") + fr(p))
    if q != 0:
        parts.append(("-" if q < 0 else "+") + (fr(q) if abs(q) != 1 else "") + r"\sqrt2")
    return parts


def mono_tex(m):
    out = ""
    for base, e in sorted(m.as_powers_dict().items(), key=lambda t: str(t[0])):
        n = names[str(base)]
        out += n + ("^{%d}" % e if e != 1 else "") if not n.startswith(r"\Li") or e == 1 else n + "^{%d}" % e
    return out


terms = sp.Poly(M3, *syms).terms()
tex = []
for mon, c in sorted(terms, key=lambda t: (-sum(t[0]), t[0]), reverse=False):
    m = sp.Mul(*[s ** e for s, e in zip(syms, mon)])
    parts = coef_tex(c)
    mt = mono_tex(m) if m != 1 else ""
    if len(parts) == 1:
        tex.append(parts[0] + mt)
    else:
        first = parts[0].lstrip("+")
        tex.append(f"+({first}{parts[1]}){mt}")
lines.append("LaTeX of M3 (grouped by monomial):")
lines.append(" ".join(tex))
for l in lines:
    print(l)
with open(os.path.join(OUTD, "m3_assemble.out"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
pickle.dump({"M3": M3, "li4s": {str(z): str(s) for z, s in li4s.items()}}, open(os.path.join(OUTD, "m3_sympy.pkl"), "wb"))
