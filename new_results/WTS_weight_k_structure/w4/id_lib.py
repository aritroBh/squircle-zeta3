"""Identification helpers.  PRECISION: call set_dps(n) AFTER importing this module (python-flint / mpmath imports are harmless,
but we never rely on import-time precision)."""
import json, os, sys
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
from w4_common import *
import flint

# ---------- exact coefficients of the 8 constants (from id_exact.py) ----------
EX = json.load(open(os.path.join(RES, "id_exact_coeffs.json")))
def _fr(s): return Fr(s)
# keys in EX look like "((2, 1), (Fraction(-1, 1), Fraction(-1, 1)))" -> parse pole (p0,p1)
import re
POLE = {(-1, -1): "-a", (0, 0): "0", (1, -1): "-b", (1, 1): "a"}
def parse_key(s):
    nums = [int(x) for x in re.findall(r"-?\d+", s)]
    i, j, _, _, _, _, p0, _, _, p1 = None, None, None, None, None, None, None, None, None, None
    # layout: (i, j), (Fraction(p0, 1), Fraction(p1, 1))  -> numbers: i, j, p0, 1, p1, 1
    i, j, p0, _, p1, _ = nums
    return f"{i}{j}|{POLE[(p0, p1)]}"
COEF = {}   # name -> {"c0":(q0,q1), "c1":..., "c2":..., "k3":...} each (rational, rational*sqrt2)
for s, d in EX.items():
    COEF[parse_key(s)] = {m: (_fr(d[m][0]), _fr(d[m][1])) for m in ("c0", "c1", "c2", "k3")}

def load_I(dps):
    d = json.load(open(os.path.join(RES, f"id_const_dps{dps}.json")))
    return {k: mp.mpf(v) for k, v in d["I"].items()}

def combo(I, which):
    """sum_n coef_n[which] * I_n as mp number (coef = q0 + q1 sqrt2)."""
    r2 = mp.sqrt(2); tot = mp.mpf(0)
    for n, c in COEF.items():
        q0, q1 = c[which]
        tot += (frac_mp(q0) + frac_mp(q1) * r2) * I[n]
    return tot

# ---------- atoms ----------
def atoms(dps_unused=None):
    r2 = mp.sqrt(2); b = r2 - 1; L = mp.log(1 + r2); lam = mp.log(2); pi2 = mp.pi**2; z3 = mp.zeta(3)
    li = mp.polylog
    Li2b = li(2, b); Li3b = li(3, b); Li3mb = li(3, -b); Li3h = li(3, 1/r2)
    A = {
     "lam^4": lam**4, "lam^3 L": lam**3*L, "lam^2 L^2": lam**2*L**2, "lam L^3": lam*L**3, "L^4": L**4,
     "pi2 lam^2": pi2*lam**2, "pi2 lam L": pi2*lam*L, "pi2 L^2": pi2*L**2,
     "Li2b lam^2": Li2b*lam**2, "Li2b lam L": Li2b*lam*L, "Li2b L^2": Li2b*L**2,
     "pi4": pi2**2, "pi2 Li2b": pi2*Li2b, "Li2b^2": Li2b**2,
     "z3 lam": z3*lam, "z3 L": z3*L,
     "Li3b lam": Li3b*lam, "Li3b L": Li3b*L, "Li3mb lam": Li3mb*lam, "Li3mb L": Li3mb*L,
     "Li3h lam": Li3h*lam, "Li3h L": Li3h*L,
     "Li4(1/2)": li(4, mp.mpf(1)/2), "Li4(b)": li(4, b), "Li4(-b)": li(4, -b), "Li4(1/sqrt2)": li(4, 1/r2),
     "Li4(-b^2)": li(4, -b**2), "Li4(2-sqrt2)": li(4, 2 - r2), "Li4(1-1/sqrt2)": li(4, 1 - 1/r2),
    }
    return A

def with_sqrt2(A):
    r2 = mp.sqrt(2); out = {}
    for k, v in A.items():
        out[k] = v; out[k + " *s2"] = v * r2
    return out

# ---------- LLL relation finder ----------
def find_relation(vals, D, delta=0.99, gram="exact"):
    """vals: list of mp numbers (first = target usually).  Returns list of candidate (coeffs, residual) from LLL rows,
    sorted by max |coeff|.  D = digits used for scaling (<= mp.mp.dps - 3)."""
    n = len(vals); sc = mp.mpf(10)**D
    rows = []
    for i, x in enumerate(vals):
        r = [0]*n; r[i] = 1; r.append(int(mp.nint(sc * x)))
        rows.append(r)
    M = flint.fmpz_mat(rows)
    R = M.lll(delta=delta, rep="zbasis", gram=gram)
    res = []
    for r in range(R.nrows()):
        co = [int(R[r, c]) for c in range(n)]
        s = sum(mp.mpf(c) * x for c, x in zip(co, vals))
        res.append((co, s))
    res.sort(key=lambda t: (abs(t[1]) > mp.mpf(10)**(-(D - 6)), max(abs(c) for c in t[0])))
    return res
