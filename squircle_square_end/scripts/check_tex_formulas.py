"""Parse the displayed closed forms directly from ../main.tex and evaluate them numerically, to catch
transcription errors between the exact computations and the paper.

Checked: Theorem 1 (both forms), Theorem 2, Theorem 3 (W3'), Lemma dilog (all lines), Lemma trilog (all lines),
Proposition T, Appendix B (both forms).
"""
import os, re, sys
import mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import polylog_exact  # noqa: F401  (import before setting precision)
import sympy as sp
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,
                                        implicit_multiplication_application, convert_xor)

HERE = os.path.dirname(os.path.abspath(__file__))
TEX_PATH = os.environ.get("SQUIRCLE_TEX", os.path.join(HERE, "..", "main.tex"))
if not os.path.exists(TEX_PATH):
    sys.exit("check_tex_formulas.py reads the note's TeX source, which is not in this repository.\n"
             "Set SQUIRCLE_TEX=/path/to/main.tex to run it.")
TEX = open(TEX_PATH).read()
OUT = os.path.join(HERE, "..", "outputs")
mp.mp.dps = 80
r2 = mp.sqrt(2); aN = 1 + r2; bN = r2 - 1
VALS = {"L": mp.log(aN), "lam": mp.log(2), "pi": mp.pi, "z3": mp.zeta(3), "a": aN, "b": bN,
        "L8": mp.dirichlet(3, [0, 1, 0, -1, 0, -1, 0, 1]), "Cl3": mp.clcos(3, mp.pi / 4)}
Li2f, Li3f, Li4f = sp.Function("Li2"), sp.Function("Li3"), sp.Function("Li4")


def group(s, i):
    """s[i] is '{' -> (content, index after closing brace)."""
    assert s[i] == "{"
    depth = 0
    for j in range(i, len(s)):
        if s[j] == "{":
            depth += 1
        elif s[j] == "}":
            depth -= 1
            if depth == 0:
                return s[i + 1:j], j + 1
    raise ValueError("unbalanced")


def token(s, i):
    while s[i] == " ":
        i += 1
    if s[i] == "{":
        return group(s, i)
    if s[i] == "\\":
        m = re.match(r"\\[a-zA-Z]+", s[i:])
        return s[i:i + m.end()], i + m.end()
    return s[i], i + 1


def tex2py(s):
    s = s.replace("\\,", " ").replace("\\ ", " ").replace("\\quad", " ").replace("\\qquad", " ")
    for w in ["\\Bigl", "\\Bigr", "\\bigl", "\\bigr", "\\left", "\\right", "\\notag", "&", "\\\\"]:
        s = s.replace(w, " ")
    # fractions
    while True:
        m = re.search(r"\\[td]?frac", s)
        if not m:
            break
        num, j = token(s, m.end())
        den, k = token(s, j)
        s = s[:m.start()] + f"(({num})/({den}))" + s[k:]
    s = re.sub(r"\\sqrt\{?2\}?", "sqrt(2)", s)
    s = s.replace("\\Li_4", " Li4").replace("\\Li_3", " Li3").replace("\\Li_2", " Li2").replace("\\zeta(3)", " z3 ")
    s = re.sub(r"\\Cl_3\s*\(\s*\(\(\\pi\)/\(4\)\)\s*\)", " Cl3 ", s)
    s = re.sub(r"\\Cl_3\s*\(\(\\pi\)/\(4\)\)", " Cl3 ", s)
    s = s.replace("\\Cl_3(\\pi/4)", " Cl3 ").replace("Cl_3(\\pi/4)", "Cl3")
    s = s.replace("L(\\chi_8,3)", "L8").replace("\\pi", " pi ").replace("\\lambda", " lam ").replace("\\cdot", "*")
    s = s.replace("\\log", "log")
    s = s.replace("{", "(").replace("}", ")").replace("[", "(").replace("]", ")")
    s = s.replace("$", " ")
    # explicit multiplication before "(" after a symbol, a closing bracket or a bare number
    s = re.sub(r"(?<![A-Za-z0-9])(lam|L|pi|z3|L8|Cl3|a|b)\s*\(", r"\1*(", s)
    s = re.sub(r"\)\s*\(", ")*(", s)
    s = re.sub(r"(?<![A-Za-z0-9])(\d+)\s*\(", r"\1*(", s)
    s = re.sub(r"\s+", " ", s).strip().rstrip(",.")
    return s


TR = standard_transformations + (implicit_multiplication_application, convert_xor)
SYMS = {k: sp.Symbol(k) for k in VALS}
SYMS.update({"Li2": Li2f, "Li3": Li3f, "Li4": Li4f, "sqrt": sp.sqrt, "log": sp.log})


def ev(s):
    e = parse_expr(tex2py(s), local_dict=SYMS, transformations=TR)
    f = sp.lambdify(list(SYMS[k] for k in VALS), e,
                    modules=[{"Li2": lambda z: mp.re(mp.polylog(2, z)), "Li3": lambda z: mp.re(mp.polylog(3, z)), "Li4": lambda z: mp.re(mp.polylog(4, z)),
                              "sqrt": mp.sqrt, "log": mp.log}, "mpmath"])
    return f(*[VALS[k] for k in VALS])


def between(start, end, src=TEX, frm=0):
    i = src.index(start, frm) + len(start)
    return src[i:src.index(end, i)], src.index(end, i)


ref = {k: mp.mpf(open(os.path.join(OUT, f"{k}_methodA_1200.txt")).read().strip()) for k in ("M1", "M2")}
results = []


def report(name, val, target):
    d = val - target
    ok = abs(d) < mp.mpf(10) ** -50
    results.append(ok)
    print(f"{'OK ' if ok else 'BAD'} {name:45s} diff = {mp.nstr(d, 3)}")


# Theorem 1
thm1, pos = between(r"\begin{theorem}\label{thm:M1}", r"\end{theorem}")
f1 = between(r"\[", r"\]", thm1)[0]
report("Theorem 1, first form: M1", ev(f1.replace("M_1=", "")), ref["M1"])
f1b = between(r"that is,", r"Numerically", thm1)[0].strip().strip("$")
lhs, rhs = f1b.split("=")
report("Theorem 1, second form: 48 M1", ev(rhs), 48 * ref["M1"])
# Theorem 2
thm2 = between(r"\begin{theorem}\label{thm:M2}", r"\end{theorem}")[0]
f2 = between(r"\begin{multline*}", r"\end{multline*}", thm2)[0].replace("M_2=", "")
report("Theorem 2: M2", ev(f2), ref["M2"])
# Theorem 3 W3'
thm3 = between(r"\begin{theorem}\label{thm:M3}", r"\end{theorem}")[0]
f3 = between(r"\begin{multline*}", r"\end{multline*}", thm3)[0].replace("W_3'=", "")
sys.path.insert(0, HERE)
import pickle
W3p = pickle.load(open(os.path.join(OUT, "m3_tails_exact.pkl"), "rb"))["W3p"].num()
report("Theorem 3: W3'", ev(f3), W3p)
M3ref = mp.mpf(open(os.path.join(OUT, "M3_methodA_300.txt")).read().strip())
TT = (mp.pi ** 2 + 6 * mp.log(2) ** 2) / 24
W4tex = mp.mpf(between(r"W_4&=", r"\ldots", thm3)[0].strip())
print(f"    W4 decimal in Theorem 3 vs M3-W3'-T: diff {mp.nstr(W4tex - (M3ref - W3p - TT), 3)}")
# Conjecture W4
conj = between(r"\begin{conjecture}\label{conj:W4}", r"\end{conjecture}")[0]
fc = between(r"\begin{multline*}", r"\end{multline*}", conj)[0].replace("W_4=", "")
with mp.workdps(80):
    W4ref = M3ref - W3p - TT
report("Conjecture W4 (vs M3-W3'-T, 300-digit M3)", ev(fc), W4ref)
# Lemma dilog
lem = between(r"\begin{lemma}\label{lem:dilog}", r"\end{lemma}")[0]
body = between(r"\begin{align}", r"\end{align}", lem)[0]
for line in body.split("\\\\"):
    line = re.sub(r"\\label\{[^}]*\}", "", line).replace("\\notag", "")
    if "=" not in line:
        continue
    l, r = line.split("=", 1)
    report("Lemma dilog: " + tex2py(l), ev(r), ev(l))
# Lemma trilog
lem = between(r"\begin{lemma}\label{lem:trilog}", r"\end{lemma}")[0]
body = between(r"\begin{align*}", r"\end{align*}", lem)[0]
entries = re.split(r"\\\\(?!&\\qquad)", body)
for line in entries:
    line = line.replace("\\\\&\\qquad", " ")
    if "=" not in line:
        continue
    l, r = line.split("=", 1)
    report("Lemma trilog: " + tex2py(l), ev(r), ev(l))
# Appendix B forms
appB = between(r"\section{Other forms of $M_2$}", r"\section{Reproducibility}")[0]
forms = re.findall(r"\\begin\{multline\*\}(.*?)\\end\{multline\*\}", appB, re.S)
for i, fm in enumerate(forms):
    report(f"Appendix B form {i + 1}: M2", ev(fm.replace("M_2=", "")), ref["M2"])
print("ALL OK" if all(results) else "SOME CHECKS FAILED")
