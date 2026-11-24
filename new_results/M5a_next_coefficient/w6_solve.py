"""Solve exactly for the depth-two correction N6 at the symbol level.

Target: W6 pairing vector T = (33: 449 sqrt2/4, 42a: -30, 42b: 15 sqrt2, 42c: -9).
Each atom has a rational pairing vector; coefficients live in Q(sqrt2).
Scan 4-subsets of a shortlist; solve the 4x4 QS2 system exactly; score by size.
"""
import itertools
import json
import sys
from fractions import Fraction as F
from pathlib import Path
import flint
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent/'A3_W4_proof'))
from w4_symbol import Q, R, A, B, ONE as ONEQ, ban_float
import w6_symbol as W

ban_float()

def vec(atom):
    (m, n), x, y = atom
    s = W.depth2_symbol(m, n, x, y)
    return [W.pair33(s)] + W.pair42(s)

def qs2(x):
    return Q(x)

SHORT = [
    ((4, 2), B, ONEQ), ((3, 3), B, ONEQ),
    ((4, 2), Q(F(0), F(1, 2)), ONEQ), ((3, 3), Q(F(0), F(1, 2)), ONEQ),
    ((4, 2), Q(2), B), ((3, 3), B, Q(2)), ((3, 3), Q(2), B),
    ((4, 2), B, Q(F(0), F(1, 2))), ((3, 3), B, Q(F(0), F(1, 2))),
    ((4, 2), Q(F(0), F(1, 2)), B), ((3, 3), Q(F(0), F(1, 2)), B),
    ((4, 2), Q(2), Q(F(0), F(1, 2))), ((3, 3), Q(2), Q(F(0), F(1, 2))),
    ((4, 2), ONEQ - Q(F(0), F(1, 2)), ONEQ), ((3, 3), ONEQ - Q(F(0), F(1, 2)), ONEQ),
    ((4, 2), Q(2) - R, ONEQ), ((3, 3), Q(2) - R, ONEQ),
    ((3, 3), -B, ONEQ), ((4, 2), -B, ONEQ),
    ((3, 3), Q(F(0), -F(1, 2)), ONEQ), ((4, 2), Q(F(0), -F(1, 2)), ONEQ),
]
NAMES = {atom: f"Li_{atom[0][0]}{atom[0][1]}({atom[1]},{atom[2]})" for atom in SHORT}

# target split into rational / sqrt2 parts
T_rat = [F(0), F(-30), F(0), F(-9)]
T_sqr = [F(449, 4), F(0), F(15), F(0)]

atomvecs = {}
for atom in SHORT:
    v = vec(atom)
    assert all(t.q == 0 for t in v), (atom, v)   # rational vectors only
    atomvecs[atom] = [F(str(t.p)) for t in v]

def solve_subset(sub):
    M = flint.fmpq_mat([[str(atomvecs[a][r]) for a in sub] for r in range(4)])
    if M.det() == 0:
        return None
    Mi = M.inv()
    pr = Mi * flint.fmpq_mat([[str(t)] for t in T_rat])
    qr = Mi * flint.fmpq_mat([[str(t)] for t in T_sqr])
    coeff = []
    size = 1
    for i in range(4):
        p = F(str(pr[i, 0])); q = F(str(qr[i, 0]))
        coeff.append((p, q))
        size = max(size, abs(p.numerator) + p.denominator, abs(q.numerator) + q.denominator)
    return coeff, size

best = []
for sub in itertools.combinations(SHORT, 4):
    r = solve_subset(sub)
    if r is None:
        continue
    coeff, size = r
    best.append((size, sub, coeff))
best.sort(key=lambda t: t[0])
print(f"solvable 4-subsets: {len(best)}")
for size, sub, coeff in best[:6]:
    print(f"size {size}:")
    for a, (p, q) in zip(sub, coeff):
        print(f"   ({p} + {q} sqrt2) * {NAMES[a]}   vec={atomvecs[a]}")

# take the best, verify against the full pairing computation
size, sub, coeff = best[0]
lhs33, lhs42 = Q(0), [Q(0)] * 3
for a, (p, q) in zip(sub, coeff):
    c = Q(p, q)
    (m, n), x, y = a
    s = W.depth2_symbol(m, n, x, y)
    lhs33 += c * W.pair33(s)
    for i, v in enumerate(W.pair42(s)):
        lhs42[i] += c * v
S6 = W.symbol(W.w6_expression())
t33 = W.pair33(S6); t42 = W.pair42(S6)
print("\nchosen solution check (exact):")
print("  N6 (3,3):", lhs33, " target:", t33, " equal:", lhs33 == t33)
print("  N6 (4,2):", [str(v) for v in lhs42], " target:", [str(v) for v in t42],
      " equal:", all(lhs42[i] == t42[i] for i in range(3)))
out = {"atoms": [], "target33": str(t33), "target42": [str(v) for v in t42]}
for a, (p, q) in zip(sub, coeff):
    out["atoms"].append({"name": NAMES[a], "coeff_p": str(p), "coeff_q": str(q),
                         "m": a[0][0], "n": a[0][1], "x": str(a[1]), "y": str(a[2])})
(HERE/'n6_solution.json').write_text(json.dumps(out, indent=1) + "\n")
print("saved n6_solution.json")
