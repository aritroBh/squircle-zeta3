import os
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
src = '''"""oct5d3_basis.py — M4f atom basis for the octagon W5 closure.

Conventions (M4e/M4f): lam = log 2 (= letter X/0), L = log(1+sqrt2) (= letter Y/1).
QS2(p,q) = p + q*sqrt2. Symbols via OCT5D2 symengine (validated w4_symbol).
Li32(x,y) value = sum_{n>m>0} x^n y^m/(n^3 m^2); its symbol is g_symbol(0,(0,0,1/x,0,1/(x*y)),1)
(OCT5D2 convention, validated: C's (3,2) pair = the obstruction).
"""
import sys
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
sys.path.insert(0, "<project>/research/08_theory/OCT5D2_octagon_depth_test")
from fractions import Fraction as F
from qsqrt2 import QS2 as Q
import symengine as S
from symengine import add, scale, shuffle, const_pow, li_symbol, gsym, obstruction

ZERO, ONE = Q(0), Q(1)
R2 = Q(0, 1)
B = Q(-1, 1)                      # b = sqrt2 - 1
H = R2/2                          # 1/sqrt2

# ---- the 14 Li5 arguments (7 W4 + 7 M4f extras) ----
LI5_ARGS = [
    ("b", B), ("-b", -B), ("-b^2", -B*B), ("1/sqrt2", H), ("1/2", Q(F(1,2))),
    ("1-1/sqrt2", ONE-H), ("2-sqrt2", 2-R2),
    ("2b", 2*B), ("(2-sqrt2)/4", (2-R2)/4), ("(2+sqrt2)/4", (2+R2)/4),
    ("-b/2", -B/2), ("(4-3sqrt2)/8", Q(F(1,2), F(-3,8))), ("b^4", B**4), ("1-b^4", ONE-B**4),
]
# ---- the 7 Li4 arguments (M4f) ----
LI4_ARGS = [("b", B), ("-b", -B), ("b^2", B*B), ("1/sqrt2", H), ("1/2", Q(F(1,2))),
            ("1-1/sqrt2", ONE-H), ("2-sqrt2", 2-R2)]
# ---- the 3 Li3 arguments ----
LI3_ARGS = [("b", B), ("-b", -B), ("1/sqrt2", H)]

def atom_list():
    """Return list of (label, kind, symbol_dict). kind: 'vis' (nonzero symbol) or 'inv' (zero symbol)."""
    atoms = []
    # 14 Li5 (visible)
    for name, z in LI5_ARGS:
        atoms.append((f"Li5({name})", 'vis', li_symbol(5, z)))
    # pure log^5: lam^i L^(5-i), i=5..0  (visible)
    for i in range(5, -1, -1):
        atoms.append((f"lam^{i} L^{5-i}" if i and 5-i else (f"lam^5" if i==5 else "L^5"),
                      'vis', shuffle(const_pow(0, i), const_pow(1, 5-i))))
    # lam^i L^j pi2, i+j=3 (invisible)
    for i in range(3, -1, -1):
        atoms.append((f"lam^{i} L^{3-i} pi2", 'inv', {}))
    # lam^i L^j Li2(b), i+j=3 (visible)
    s2b = li_symbol(2, B)
    for i in range(3, -1, -1):
        atoms.append((f"lam^{i} L^{3-i} Li2(b)", 'vis',
                      shuffle(s2b, shuffle(const_pow(0, i), const_pow(1, 3-i)))))
    # lam pi2^2, L pi2^2 (invisible)
    atoms.append(("lam pi2^2", 'inv', {})); atoms.append(("L pi2^2", 'inv', {}))
    # lam pi2 Li2(b), L pi2 Li2(b) (invisible)
    atoms.append(("lam pi2 Li2(b)", 'inv', {})); atoms.append(("L pi2 Li2(b)", 'inv', {}))
    # lam Li2(b)^2, L Li2(b)^2 (visible)
    s22 = shuffle(s2b, s2b)
    atoms.append(("lam Li2(b)^2", 'vis', shuffle(s22, const_pow(0, 1))))
    atoms.append(("L Li2(b)^2", 'vis', shuffle(s22, const_pow(1, 1))))
    # lam^i L^j z3, i+j=2 (invisible)
    for i in range(2, -1, -1):
        atoms.append((f"lam^{i} L^{2-i} z3", 'inv', {}))
    # lam^i L^j Li3(z), i+j=2 (visible)
    for i in range(2, -1, -1):
        for name, z in LI3_ARGS:
            atoms.append((f"lam^{i} L^{2-i} Li3({name})", 'vis',
                          shuffle(li_symbol(3, z), shuffle(const_pow(0, i), const_pow(1, 2-i)))))
    # pi2 z3 (invisible)
    atoms.append(("pi2 z3", 'inv', {}))
    # pi2 Li3(z) (invisible)
    for name, z in LI3_ARGS:
        atoms.append((f"pi2 Li3({name})", 'inv', {}))
    # Li2(b) z3 (invisible)
    atoms.append(("Li2(b) z3", 'inv', {}))
    # Li2(b) Li3(z) (visible)
    for name, z in LI3_ARGS:
        atoms.append((f"Li2(b) Li3({name})", 'vis', shuffle(s2b, li_symbol(3, z))))
    # lam Li4(z), L Li4(z) (visible)
    for name, z in LI4_ARGS:
        s4 = li_symbol(4, z)
        atoms.append((f"lam Li4({name})", 'vis', shuffle(s4, const_pow(0, 1))))
        atoms.append((f"L Li4({name})", 'vis', shuffle(s4, const_pow(1, 1))))
    # z5, Lchi8_5 (invisible)
    atoms.append(("z5", 'inv', {})); atoms.append(("Lchi8_5", 'inv', {}))
    return atoms

# ---- the depth-two correction C (fixed coefficients) ----
C_TERMS = [("Li32(b,1)", (B, ONE), Q(-3600, 2160)),
           ("Li32(1/sqrt2,1)", (H, ONE), Q(11880, -8280))]

def w32(x, y):
    return (ZERO, ZERO, ONE/x, ZERO, ONE/(x*y))

def li32_symbol(x, y):
    return gsym(w32(x, y))

def C_symbol():
    out = {}
    for name, (x, y), cf in C_TERMS:
        out = add(out, scale(li32_symbol(x, y), cf))
    return out

if __name__ == "__main__":
    import time
    t0 = time.time()
    atoms = atom_list()
    vis = [a for a in atoms if a[1] == 'vis']
    inv = [a for a in atoms if a[1] == 'inv']
    print("atoms:", len(atoms), "visible:", len(vis), "invisible:", len(inv),
          "built in %.1fs" % (time.time()-t0))
    assert len(atoms) == 70 and len(vis) == 52 and len(inv) == 18
    cs = C_symbol()
    print("C pair (should be the obstruction (-1485+1035sqrt2, 900-540sqrt2)):",
          obstruction(cs))
'''
open("oct5d3_basis.py","w").write(src)
print("written")
import subprocess
r = subprocess.run([sys.executable, "oct5d3_basis.py"], capture_output=True, text=True)
print(r.stdout); print(r.stderr[-2000:] if r.returncode else "OK")
