"""Per-value table for Lemma (trilog): for each Li3 value occurring in W3 (M2) and W3' (M3), express it in the
basis Li3(b), Li3(-b), Li3(1/sqrt2), zeta(3) + products, and record whether the one-variable relations
(I3), (T3), (D3), (V3) suffice or the Spence-Kummer instances (K) are needed.
"""
import os, sys, pickle
import mpmath as mp
from fractions import Fraction as Fr
from qsqrt2 import QS2, SQ2, A, B
import polylog_exact as PE
from polylog_exact import P, LI3, ZERO, ONE
import kummer

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 60
atoms = [QS2(-2, -2), -A, QS2(-1), QS2(0, Fr(-1, 2)), -B / 2, -B * B, (QS2(2) - SQ2) / 4, ONE - ONE / SQ2,
         QS2(Fr(1, 2)), QS2(2) - SQ2, B * 2, (QS2(2) + SQ2) / 4]
basis = [B, -B, QS2(0, Fr(1, 2))]
KREL = [(f"Kummer (K) at (x,y)=({x}, {y})", kummer.kummer_relation(x, y))
        for (x, y) in [(B, QS2(0, Fr(1, 2))), (QS2(2) - SQ2, QS2(0, Fr(1, 2))), (QS2(0, Fr(1, 2)), 2 * B), (QS2(-1), B * B)]]
seed = set(atoms) | set(basis)
CLASSICAL = PE.li3_closure(seed, 6, ())
for nm, rel in CLASSICAL + KREL:
    PE._check(nm, rel, P())


def solve(expr, rels):
    allatoms = sorted({a for _, r in rels for a in r.atoms('Li3')} | expr.atoms('Li3'), key=float)
    nonb = [a for a in allatoms if a not in basis]
    idx = {a: i for i, a in enumerate(nonb)}

    def vec_of(p):
        v = [ZERO] * len(nonb)
        for k, c in p.d.items():
            if len(k) == 1 and not isinstance(k[0], str) and k[0][0] == 'Li3' and k[0][1] in idx:
                v[idx[k[0][1]]] = v[idx[k[0][1]]] + c
        return v
    cols = [vec_of(r) for _, r in rels]
    tgt = vec_of(expr)
    n, m = len(nonb), len(cols)
    Mx = [[cols[j][i] for j in range(m)] + [tgt[i]] for i in range(n)]
    piv, row = [], 0
    for col in range(m):
        pr = next((r for r in range(row, n) if Mx[r][col] != ZERO), None)
        if pr is None:
            continue
        Mx[row], Mx[pr] = Mx[pr], Mx[row]
        inv = ONE / Mx[row][col]
        Mx[row] = [x * inv for x in Mx[row]]
        for r in range(n):
            if r != row and Mx[r][col] != ZERO:
                f = Mx[r][col]
                Mx[r] = [x - f * y for x, y in zip(Mx[r], Mx[row])]
        piv.append(col); row += 1
        if row == n:
            break
    if any(Mx[i][m] != ZERO for i in range(row, n)):
        return None
    red = expr
    for r_, col in enumerate(piv):
        if Mx[r_][m] != ZERO:
            red = red - rels[col][1] * Mx[r_][m]
    assert not [a for a in red.atoms('Li3') if a not in basis]
    return red


lines = []
for z in atoms:
    red = solve(LI3(z), CLASSICAL)
    kind = "classical (I3,T3,D3,V3)"
    if red is None:
        red = solve(LI3(z), CLASSICAL + KREL)
        kind = "needs Spence-Kummer (K)"
    with mp.workdps(60):
        err = LI3(z).num() - red.num()
    s = f"Li3({z}) = {red}    [{kind}; numeric check {mp.nstr(err, 3)}]"
    print(s, flush=True)
    lines.append(s)
with open(os.path.join(HERE, "..", "outputs", "li3_table.out"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
