
# step_06b.py -- delta22 of Omega as a FUNCTION of c, per exponent group, vs paper Prop. M3angle
import sys; sys.path.insert(0, '.')
from cor6x_core import *
import sympy as sp, json, re
c = sp.Symbol('c')
data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
TERMS = [(parse_term(s)[0], parse_term(s)[1], sp.Rational(cf)) for s, cf in data if len(parse_term(s)[1]) == 4]

# letters modulo constants: map each letter expr to a representative in {c, 1-c, 1+c} or None
def modconst(S):
    out = {}
    for mono, v in S.items():
        new = []
        ok = True
        for t in mono:
            tt = sp.expand(t)
            if tt == c:
                new.append(c)
            elif sp.simplify(tt - (1 - c)) == 0 or sp.simplify(tt - (c - 1)) == 0:
                new.append(1 - c)
            elif sp.simplify(tt - (1 + c)) == 0 or sp.simplify(tt - (-1 - c)) == 0:
                new.append(1 + c)
            else:
                ok = False
                break
        if ok and len(new) == len(mono):
            out[tuple(new)] = out.get(tuple(new), 0) + v
    return {k: v for k, v in out.items() if v}

groups = {}
for e, w, cf in TERMS:
    groups.setdefault(e, []).append((w, cf))
print("exponent groups:", {e: len(v) for e, v in groups.items()})
Dtot = {}
for e in sorted(groups):
    Sg = {}
    for w, cf in groups[e]:
        Sg = sym_add(Sg, sym_scale(sym_H(w, c), cf))
    # check which letters appear before mod-const
    letters = set(t for mono in Sg for t in mono)
    print(f"e={e}: symbol terms {len(Sg)}, letters: {sorted(map(str, letters))}")
    Dg = delta22(modconst(Sg))
    print(f"   delta22 (mod constants): {Dg}")
    for k, v in Dg.items():
        Dtot[k] = Dtot.get(k, 0) + c**e * v
print()
print("delta22(Omega)(c) = sum_e c^e * group:", Dtot)
# paper: (1-1/c)([c]^[1-c])^([1-c]^[1+c]) - (1+1/c)([c]^[1+c])^([1-c]^[1+c])
