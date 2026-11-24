"""Split M_k into pure pieces of uniform weight by integrating the higher poles at v = +-i by parts.

M_k = int_b^a sum_{(i,j)} R_ij(v) l1^i l2^j dv,  l1 = log t, l2 = log(1-t)   (mk_integrands.py)
Each R_ij = (simple poles at real letters 0, -a, +-b, a) + (poles of order >= 2 at +-i) [+ simple poles at +-i].
The order>=2 part H is exact: H = G', G rational. With C = G(b),
    int l1^i l2^j H dv = - int (G - C) d(l1^i l2^j)          (boundary terms vanish, see note)
which lowers the log-degree by one. Monomials with i < j are folded by v -> 1/v (t -> 1-t):
    int l1^i l2^j R(v) dv = int l1^j l2^i R(1/v) v^-2 dv.
Iterate until only simple poles remain. Result: {(i,j): R_ij} with simple poles only;
the piece with i+j = m is a pure iterated integral of weight m+1.
"""
import sympy as sp
from mk_integrands import v_integrand, v, S2

a_, b_ = 1 + S2, S2 - 1
I = sp.I
DL1 = 1 / (v - b_) + 1 / (v + a_) - 1 / v        # d log t / dv
DL2 = 1 / (v - a_) + 1 / (v + b_) - 1 / v        # d log(1-t) / dv


def split(R):
    """R -> (simple-pole part, higher-order part at +-i, polynomial part)."""
    pf = sp.apart(sp.together(R), v, extension=[S2, I])
    simple, higher, poly = 0, 0, 0
    for term in sp.Add.make_args(sp.expand(pf)):
        num, den = sp.fraction(sp.together(term))
        if not den.has(v):
            poly += term
            continue
        dp = sp.Poly(den, v)
        roots = sp.roots(dp)
        assert len(roots) == 1, term
        (root, mult), = roots.items()
        if mult >= 2:
            assert sp.simplify(root ** 2 + 1) == 0, f"higher pole at {root}"
            higher += term
        else:
            simple += term
    return sp.simplify(simple), higher, sp.simplify(poly)


def fold(key, R):
    i, j = key
    if i >= j:
        return key, R
    return (j, i), sp.together(R.subs(v, 1 / v) / v ** 2)


def reduce_k(k, verbose=True):
    work = {}
    for key, R in v_integrand(k).items():
        kk, RR = fold(key, R)
        work[kk] = sp.together(work.get(kk, 0) + RR)
    final = {}
    for deg in range(k, -1, -1):
        for key in sorted([kk for kk in work if sum(kk) == deg], reverse=True):
            R = work.pop(key)
            if deg == 0:
                final[key] = sp.together(sp.radsimp(R))
                continue
            simple, higher, poly = split(R)
            assert poly == 0, f"polynomial part at {key}: {poly}"
            if simple != 0:
                final[key] = sp.together(final.get(key, 0) + simple)
            if higher == 0:
                continue
            G = sp.integrate(higher, v)
            G = sp.together(G)
            assert sp.simplify(sp.diff(G, v) - higher) == 0
            C = sp.radsimp(sp.simplify(G.subs(v, b_)))
            Gm = G - C
            i, j = key
            if i >= 1:
                kk, RR = fold((i - 1, j), -i * Gm * DL1)
                work[kk] = sp.together(work.get(kk, 0) + RR)
            if j >= 1:
                kk, RR = fold((i, j - 1), -j * Gm * DL2)
                work[kk] = sp.together(work.get(kk, 0) + RR)
    out = {}
    for key, R in final.items():
        R = sp.radsimp(sp.together(R))
        if key == (0, 0):
            out[key] = sp.simplify(R)
        else:
            pf = sp.apart(R, v, extension=[S2, I])
            out[key] = pf
    if verbose:
        print(f"===== k = {k}: pure pieces (weight = i+j+1)")
        for key in sorted(out, key=lambda kk: (-sum(kk), kk)):
            print(f"  weight {sum(key)+1}  l1^{key[0]} l2^{key[1]}:")
            for term in sp.Add.make_args(out[key]):
                print("        ", sp.simplify(term))
    return out


if __name__ == "__main__":
    import sys
    for k in [int(x) for x in sys.argv[1:]] or [2, 3]:
        reduce_k(k)
