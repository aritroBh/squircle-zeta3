"""Independent known values (closed forms from the paper / repo).  Call after set_dps()."""
import re, json, os
from w4_common import *
from fractions import Fraction as F

def consts():
    r2 = mp.sqrt(2); b = r2 - 1; L = mp.log(1 + r2); lam = mp.log(2)
    return r2, b, L, lam

def W2_exact():
    r2, b, L, lam = consts()
    return (mp.pi**2 + 6 * lam**2) / 24

def S_const():
    r2, b, L, lam = consts()
    return mp.zeta(3) / 8 + mp.pi**2 * lam / 24 + lam**3 / 12

def W3p():
    """W3' of paper/squircle_note/main.tex thm:M3 (weight-3 part of M_3), transcribed from the multline."""
    r2, b, L, lam = consts()
    li = mp.polylog
    return ((7*r2/32 - mp.mpf(3)/16) * mp.zeta(3) + (1 + r2) * li(3, b) + (1 - r2) * li(3, -b)
            - 2*r2 * li(3, 1/r2) + (2*L - r2*lam) * li(2, b)
            + (r2/12 - mp.mpf(1)/8) * mp.pi**2 * L - (1 + r2)/24 * mp.pi**2 * lam + (1 + r2)/24 * lam**3
            + r2/4 * lam**2 * L - (1 + r2)/2 * lam * L**2 + (5 - r2)/6 * L**3)

def W3_known(k):
    return W3p() - (k - 3) * S_const()

def W4_thm():
    """W_4 of Theorem thm:W4 (main.tex), transcribed term by term from the multline."""
    r2, b, L, lam = consts()
    li = mp.polylog; z3 = mp.zeta(3); pi2 = mp.pi**2; q = mp.mpf
    t = ((1 - 573*r2) * li(4, b) - (1 + 573*r2) * li(4, -b) + 56*r2 * li(4, -b**2)
         + 17*r2 * li(4, 1/r2) - (q(37)/4 + 17*r2/16) * li(4, q(1)/2)
         + (10 + 86*r2) * li(4, 1 - 1/r2) - (10 - 86*r2) * li(4, 2 - r2))
    t += L * ((40 + 7*r2) * li(3, b) - (40 - 7*r2) * li(3, -b) - 72 * li(3, 1/r2) + (q(63)/8 - 3*r2/16) * z3)
    t += lam * ((6 + 13*r2) * li(3, b) + (6 - 13*r2) * li(3, -b) - 47*r2/2 * li(3, 1/r2) + (329*r2/128 - q(3)/8) * z3)
    t += 30 * li(2, b)**2
    t += li(2, b) * ((30 + 7*r2/2) * L**2 - 24*lam*L - 19*r2/8 * lam**2 + (r2/24 - q(19)/4) * pi2)
    t += (q(103)/384 - 4963*r2/5760) * pi2**2
    t += pi2 * ((q(67)/24 + 185*r2/96) * L**2 + (41*r2/12 - q(35)/12) * lam*L + (q(5)/12 - 1105*r2/384) * lam**2)
    t += ((q(13)/24 - 25*r2/24) * L**4 - (q(73)/6 + 29*r2/3) * L**3 * lam + (q(63)/8 + 67*r2/16) * L**2 * lam**2
          + (q(41)/24 + 17*r2/6) * L * lam**3 + (137*r2/192 - q(61)/192) * lam**4)
    return t

def W4pp_txt():
    """numeric W4'' = weight-4 part of M_4 from M4a_decomposition/W4pp_800.txt (800 digits)."""
    p = os.path.join(os.path.dirname(HERE), "..", "M4a_decomposition", "W4pp_800.txt")
    p = os.path.normpath(p)
    s = re.sub(r"\s", "", open(p).read()).split("=")[-1]
    return mp.mpf(s)

def _field(s):
    r2 = mp.sqrt(2)
    return eval(re.sub(r"(\d+)", r"mp.mpf(\1)", s.replace("sqrt2", "RT")), {"mp": mp, "RT": r2})

def _split(name):
    out, depth, cur = [], 0, ""
    for ch in name:
        if ch == "(": depth += 1
        if ch == ")": depth -= 1
        if ch == "*" and depth == 0: out.append(cur); cur = ""
        else: cur += ch
    out.append(cur); return out

def W4pp_closed():
    """W4'' evaluated from w4pp_closed_form.json with my own atom evaluator (independent of planner_eval_R5)."""
    r2, b, L, lam = consts()
    p = os.path.normpath(os.path.join(os.path.dirname(HERE), "..", "M4a_decomposition", "w4pp_closed_form.json"))
    terms = json.load(open(p))["closed_form_terms"]
    def atom(a):
        if a == "sqrt2": return r2
        if a == "log2": return lam
        if a == "L": return L
        if a == "pi^2": return mp.pi**2
        if a == "z3": return mp.zeta(3)
        if a == "Li2(b)": return mp.polylog(2, b)
        if a == "Li3(b)": return mp.polylog(3, b)
        if a == "Li3(-b)": return mp.polylog(3, -b)
        if a == "Li3(1/sqrt2)": return mp.polylog(3, 1/r2)
        m = re.fullmatch(r"Li4\((.*)\)", a)
        assert m, a
        return mp.polylog(4, _field(m.group(1)))
    tot = mp.mpf(0)
    for c, name in terms:
        v = frac_mp(F(c))
        for f in _split(name): v *= atom(f)
        tot += v
    return tot
