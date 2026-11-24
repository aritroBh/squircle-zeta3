"""Assemble the W4'' candidate closed form from the LLL relation and verify numerically."""
import sys, json
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import mpmath as mp
from fractions import Fraction as Fr
from qsqrt2 import QS2, SQ2, A, B
import polylog_exact as PE

mp.mp.dps = 1120
rel = json.load(open("pslq2_w4pp_700.json"))["relation"]
c0 = rel["W4''"]
assert c0 == 11520
W = mp.mpf(open("W4pp_1060.txt").read().strip())

r2 = mp.sqrt(2); lam = mp.log(2); L = mp.log(1+r2); b = r2-1
consts = {"log2": lam, "L": L, "pi^2": mp.pi**2, "Li2(b)": mp.polylog(2, b),
          "z3": mp.zeta(3), "Li3(b)": mp.polylog(3, b), "Li3(-b)": mp.polylog(3, -b),
          "Li3(1/sqrt2)": mp.polylog(3, 1/r2)}
ZMAP = {"Li4(1/2)": QS2(Fr(1,2)), "Li4((-1+1*sqrt2))": B, "Li4((1-1*sqrt2))": -B,
        "Li4(1/2*sqrt2)": QS2(0, Fr(1,2)), "Li4((-3+2*sqrt2))": -B*B,
        "Li4((2-1*sqrt2))": QS2(2)-SQ2, "Li4((1-1/2*sqrt2))": QS2(1)-QS2(0, Fr(1,2))}
def val_of(name):
    if name.startswith("Li4("):
        return mp.polylog(4, PE.DE.float_exact(ZMAP[name]))
    if name.startswith("sqrt2*Li4("):
        return r2*mp.polylog(4, PE.DE.float_exact(ZMAP[name[6:]]))
    if name.startswith("sqrt2*"):
        return r2*mp.fprod(consts[p] for p in name[6:].split("*"))
    return mp.fprod(consts[p] for p in name.split("*"))

tot = mp.mpf(0)
for nm, c in rel.items():
    if nm == "W4''": continue
    tot += c*val_of(nm)
Wcf = -tot/c0
print("closed form - W4''(1060-digit quadrature) =", mp.nstr(Wcf - W, 6))

# Pretty exact form: W4'' = sum -(c/11520) name
lines = []
def fmt(name, c):
    frac = Fr(-c, c0)
    return frac, name
items = [(fmt(nm, c)[0], nm) for nm, c in rel.items() if nm != "W4''"]
json.dump({"closed_form_terms": [[str(f), nm] for f, nm in items]},
          open("w4pp_closed_form.json", "w"), indent=1)
# print grouped readable
for f, nm in items:
    print(f"  {f} * {nm}")
