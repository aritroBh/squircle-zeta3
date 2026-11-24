"""EXACT (rational-arithmetic) consistency check of the LLL relation for c_0:
   c_0 = W_4(4) - 4 c_1 - 16 c_2, with W_4(4) = W4'' from M4a_decomposition/w4pp_closed_form.json (49 terms, exact rationals),
   c_1, c_2 from the LLL relations.  The coordinates in the 58-dim Q-basis must agree EXACTLY with the LLL relation for c_0.
   Also: the same for W_4(3)-polynomial offset is NOT checked here (delta is not in the basis claim)."""
import json, os, re, sys
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))
norm = json.load(open(os.path.join(HERE, "results", "id_lll_normalised.json")))
get = lambda nm: {k: F(v) for k, v in next(iter(norm[nm].values())).items()}
c0, c1, c2 = get("c0"), get("c1"), get("c2")
terms = json.load(open(os.path.join(HERE, "..", "..", "M4a_decomposition", "w4pp_closed_form.json")))["closed_form_terms"]
def split(name):
    out, depth, cur = [], 0, ""
    for ch in name:
        if ch == "(": depth += 1
        if ch == ")": depth -= 1
        if ch == "*" and depth == 0: out.append(cur); cur = ""
        else: cur += ch
    out.append(cur); return out
LI4 = {"Li4((-1+1*sqrt2))": "Li4(b)", "Li4((1-1*sqrt2))": "Li4(-b)", "Li4(1/2*sqrt2)": "Li4(1/sqrt2)",
       "Li4((-3+2*sqrt2))": "Li4(-b^2)", "Li4((2-1*sqrt2))": "Li4(2-sqrt2)", "Li4((1-1/2*sqrt2))": "Li4(1-1/sqrt2)", "Li4(1/2)": "Li4(1/2)"}
def label(facs):
    cnt = {}
    s2 = False
    for f in facs:
        if f == "sqrt2": s2 = True; continue
        cnt[f] = cnt.get(f, 0) + 1
    def pw(n, e): return "" if e == 0 else (f"{n}^{e}" if e > 1 else n)
    lam, L = cnt.pop("log2", 0), cnt.pop("L", 0)
    pi2 = cnt.pop("pi^2", 0); li2 = cnt.pop("Li2(b)", 0); z3 = cnt.pop("z3", 0)
    li3 = {k: cnt.pop(k, 0) for k in ("Li3(b)", "Li3(-b)", "Li3(1/sqrt2)")}
    li4 = [k for k in list(cnt) if k.startswith("Li4")]
    if li4:
        assert len(cnt) == 1 and cnt[li4[0]] == 1; base = LI4[li4[0]]
    else:
        assert not cnt, cnt
        # canonical labels used in id_lib.atoms()
        parts = {
          (4,0,0,0,0,0): "lam^4", (3,1,0,0,0,0): "lam^3 L", (2,2,0,0,0,0): "lam^2 L^2", (1,3,0,0,0,0): "lam L^3", (0,4,0,0,0,0): "L^4",
        }
        key = (lam, L, pi2, li2, z3, sum(li3.values()))
        if key in parts: base = parts[key]
        elif pi2 == 1 and li2 == 0 and z3 == 0 and sum(li3.values()) == 0 and lam + L == 2:
            base = {(2,0): "pi2 lam^2", (1,1): "pi2 lam L", (0,2): "pi2 L^2"}[(lam, L)]
        elif li2 == 1 and pi2 == 0 and z3 == 0 and sum(li3.values()) == 0 and lam + L == 2:
            base = {(2,0): "Li2b lam^2", (1,1): "Li2b lam L", (0,2): "Li2b L^2"}[(lam, L)]
        elif pi2 == 2 and lam + L == 0 and li2 == 0 and z3 == 0 and sum(li3.values()) == 0: base = "pi4"
        elif pi2 == 1 and li2 == 1 and lam + L == 0 and z3 == 0: base = "pi2 Li2b"
        elif li2 == 2 and pi2 == 0 and lam + L == 0 and z3 == 0: base = "Li2b^2"
        elif z3 == 1 and lam + L == 1 and sum(li3.values()) == 0 and pi2 == 0 and li2 == 0: base = "z3 lam" if lam else "z3 L"
        elif sum(li3.values()) == 1 and lam + L == 1 and pi2 == 0 and li2 == 0 and z3 == 0:
            n = [k for k, v in li3.items() if v][0]
            nm = {"Li3(b)": "Li3b", "Li3(-b)": "Li3mb", "Li3(1/sqrt2)": "Li3h"}[n]
            base = nm + (" lam" if lam else " L")
        else: raise ValueError((facs, key))
    return base + (" *s2" if s2 else "")
W44 = {}
for c, name in terms:
    lab = label(split(name)); W44[lab] = W44.get(lab, F(0)) + F(c)
expected = {}
for lab in set(W44) | set(c1) | set(c2):
    v = W44.get(lab, F(0)) - 4 * c1.get(lab, F(0)) - 16 * c2.get(lab, F(0))
    if v: expected[lab] = v
print("W4'' atoms:", len(W44), " expected c0 atoms:", len(expected), " LLL c0 atoms:", len(c0))
print("EXACT AGREEMENT  c0(LLL) == W4'' - 4 c1(LLL) - 16 c2(LLL):", expected == c0)
if expected != c0:
    for k in sorted(set(expected) | set(c0)):
        if expected.get(k, 0) != c0.get(k, 0): print("  diff", k, expected.get(k, 0), c0.get(k, 0))
