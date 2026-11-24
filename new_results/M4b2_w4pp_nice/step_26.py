base = '<project>'
here = base + '/research/08_theory/M4b2_w4pp_nice'
code = r'''"""M4b2 step10: emit the readable Theorem-W4-style LaTeX display of R''_nice and the
M4a comparison table (markdown) from stage1_exact_cz_nice.json. Writes nice_display.md."""
import sys, json
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/research/08_theory/M4b2_w4pp_nice")
import m4b2_common as mc
from fractions import Fraction as F

d = json.load(open(str(mc.HERE/"stage1_exact_cz_nice.json")))
names = d["atom_names"]
cr = [F(x) for x in d["c_rat"]]; cs = [F(x) for x in d["c_sq"]]

LI4TEX = {"Li4((-1+1*sqrt2))": r"\Li_4(b)", "Li4((1-1*sqrt2))": r"\Li_4(-b)",
          "Li4((-3+2*sqrt2))": r"\Li_4(-b^2)", "Li4(1/2*sqrt2)": r"\Li_4\bigl(\tfrac1{\sqrt2}\bigr)",
          "Li4(1/2)": r"\Li_4\bigl(\tfrac12\bigr)", "Li4((1-1/2*sqrt2))": r"\Li_4\bigl(1-\tfrac1{\sqrt2}\bigr)",
          "Li4((2-1*sqrt2))": r"\Li_4(2-\sqrt2)"}
PRODTEX = {"log2*log2*log2*log2": r"\lambda^4", "log2*log2*log2*L": r"\lambda^3 L",
           "log2*log2*L*L": r"\lambda^2 L^2", "log2*L*L*L": r"\lambda L^3", "L*L*L*L": r"L^4",
           "pi^2*log2*log2": r"\pi^2\lambda^2", "pi^2*log2*L": r"\pi^2\lambda L",
           "pi^2*L*L": r"\pi^2 L^2", "Li2(b)*log2*log2": r"\Li_2(b)\lambda^2",
           "Li2(b)*log2*L": r"\Li_2(b)\lambda L", "Li2(b)*L*L": r"\Li_2(b)L^2",
           "pi^2*pi^2": r"\pi^4", "pi^2*Li2(b)": r"\pi^2\Li_2(b)", "Li2(b)*Li2(b)": r"\Li_2(b)^2",
           "z3*log2": r"\zeta(3)\lambda", "z3*L": r"\zeta(3)L",
           "Li3(b)*log2": r"\Li_3(b)\lambda", "Li3(b)*L": r"\Li_3(b)L",
           "Li3(-b)*log2": r"\Li_3(-b)\lambda", "Li3(-b)*L": r"\Li_3(-b)L",
           "Li3(1/sqrt2)*log2": r"\Li_3\bigl(\tfrac1{\sqrt2}\bigr)\lambda",
           "Li3(1/sqrt2)*L": r"\Li_3\bigl(\tfrac1{\sqrt2}\bigr)L"}

def fr(x):
    if x.denominator == 1: return str(x.numerator)
    return ("-" if x < 0 else "") + "\\tfrac{%d}{%d}" % (abs(x.numerator), x.denominator)

def coef_tex(p, q):
    if not q: return fr(p)
    if not p:
        if q.denominator == 1 and abs(q.numerator) == 1: return ("" if q > 0 else "-") + "\\sqrt2"
        return fr(q) + "\\sqrt2"
    sgn = "+" if q > 0 else "-"
    qa = abs(q)
    qs = ("" if qa.denominator == 1 and qa.numerator == 1 else fr(qa)) + "\\sqrt2"
    return "\\bigl(%s%s%s\\bigr)" % (fr(p), sgn, qs)

def atom_tex(name):
    return LI4TEX.get(name) or PRODTEX[name]

def term_tex(p, q, name, first=False):
    assert p or q
    neg = (p < 0) or (p == 0 and q < 0)
    if neg: p, q = -p, -q
    one = (q == 0 and p == 1)
    t = atom_tex(name) if one else coef_tex(p, q) + r"\," + atom_tex(name)
    if first and not neg: return t
    return ("- " if neg else "+ ") + t

li4 = [(k, n) for k, n in enumerate(names) if n.startswith("Li4")]
pr = [(k, n) for k, n in enumerate(names) if not n.startswith("Li4")]
lines1 = [term_tex(cr[k], cs[k], n, first=(i == 0)) for i, (k, n) in enumerate(li4)]
lines2 = [term_tex(cr[k], cs[k], n) for (k, n) in pr]
NL = chr(10)
tex = "W_4'' = " + (" \\\\" + NL).join(lines1) + " \\\\[2pt]" + NL + (" \\\\" + NL).join(lines2)

rows = ["| atom $B_k$ | $c_{\\rm rat}$ | $c_{\\rm sq}$ |", "|---|---|---|"]
for k, n in enumerate(names):
    rows.append("| %s | %s | %s |" % (n, cr[k], cs[k]))
md = (["## R''_nice in readable (Theorem W4 style) form", "",
       "With $\\lambda=\\log 2$, $L=\\log(1+\\sqrt2)$, $b=\\sqrt2-1$ (and $\\Li_s$ the classical polylogarithm):", "",
       "$$", tex, "$$", "",
       "## Coefficient table (each $c_k = c_{\\rm rat} + \\sqrt2\\,c_{\\rm sq}$; identical to M4a w4pp_closed_form.json, 58/58 components)", ""]
      + rows)
open(str(mc.HERE/"nice_display.md"), "w").write(NL.join(md))
mc.log("wrote nice_display.md")
'''
open(here + '/step10_emit_display.py', 'w').write(code)
import subprocess
PY = base + '/.venv/bin/python'
r = subprocess.run(['nice', PY, 'step10_emit_display.py'], cwd=here, capture_output=True, text=True)
print('rc=', r.returncode, r.stderr[-500:])
print(open(here + '/nice_display.md').read()[:3300])
