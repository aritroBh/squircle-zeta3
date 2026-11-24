# Copied from research/08_theory/REF_HEX5_referee/oct_weight5_fast.py (referee's code). Changes: module name, output directory $G9WORK.
# The notes' table parsed_mons.pkl is only read, to compare.
import os
"""Same computation as oct_weight5_integrand.py (leading-log = weight-5 part of the regular octagon's
order-five edge integrand), in exact Q(sqrt2) arithmetic with dict polynomials (fast).
Polynomials: dict (e0,e1,e2,e3,n) -> Q2, e_j = power of sigma_j = log T_j(t), n = power of t.
Output: rows {(a,b,c,d,e,f): {n: Q2}} in OCT5D's normalisation (P^{(5)}(1)/u = 960 * sum), saved to
own_weight5_rows.pkl, and compared with OCT5D's parsed_mons.pkl (read only)."""
import pickle, math, itertools
from fractions import Fraction as Fr
from g9_oct5_symbol import Q2, unitvec, ONE, ZERO

NV = 4
DEG = 5


def padd(p, q, c=None):
    out = dict(p)
    for k, v in q.items():
        v2 = v if c is None else v * c
        if k in out:
            s = out[k] + v2
            if s.iszero():
                del out[k]
            else:
                out[k] = s
        elif not v2.iszero():
            out[k] = v2
    return out


def pmul(p, q):
    out = {}
    for k1, v1 in p.items():
        d1 = sum(k1[:NV])
        for k2, v2 in q.items():
            if d1 + sum(k2[:NV]) > DEG:
                continue
            k = tuple(a + b for a, b in zip(k1, k2))
            s = out.get(k)
            s = v1 * v2 if s is None else s + v1 * v2
            if s.iszero():
                out.pop(k, None)
            else:
                out[k] = s
    return out


def pscale(p, c):
    return {k: v * c for k, v in p.items() if not (v * c).iszero()}


def const(c):
    return {(0,) * (NV + 1): c if isinstance(c, Q2) else Q2(c)}


def tpoly(c0, c1):
    out = {}
    if not c0.iszero():
        out[(0,) * (NV + 1)] = c0
    if not c1.iszero():
        out[(0,) * NV + (1,)] = c1
    return out


def main():
    r2 = Q2(0, 1)
    R = Q2(-1, 1)
    h = Q2(0, Fr(1, 2))  # sqrt2/2
    Vs = [(ZERO, R), (-h * R, h * R), (-R, ZERO), (-h * R, -h * R)]
    nvec = [(-v[1] / R, v[0] / R) for v in Vs]
    dot = lambda a, b: a[0] * b[0] + a[1] * b[1]
    for q in range(4):
        s = Q2(0)
        for j in range(4):
            x = dot(nvec[j], Vs[q])
            s = s + (x if float(x.a) + float(x.b) * math.sqrt(2) >= 0 else -x)
        assert s == ONE
    P0, P1 = Vs[0], Vs[1]
    D = (P1[0] - P0[0], P1[1] - P0[1])
    L2 = dot(D, D); q0 = dot(P0, D); V02 = dot(P0, P0)
    T = []
    for j in range(4):
        v0 = dot(nvec[j], P0); v1 = dot(nvec[j], P1)
        sg = 1 if float((v0 + v1).a) + float((v0 + v1).b) * math.sqrt(2) > 0 else -1
        c0 = v0 * Q2(sg); c1 = (v1 - v0) * Q2(sg)
        T.append((c0, c1))
    tot0 = sum((c[0] for c in T), Q2(0)); tot1 = sum((c[1] for c in T), Q2(0))
    assert tot0 == ONE and tot1.iszero()
    info = []
    for (c0, c1) in T:
        if c0.iszero():
            info.append(('a', c0 + c1))
        elif (c0 + c1).iszero():
            info.append(('b', c0))
        else:
            info.append(('nl', c0, c1 / c0))
    print('terms on the edge (kind, constant, kappa):', info)
    kap0 = Q2(-1, Fr(1, 2)); kap1 = Q2(-1, 1)
    # series
    unit = lambda j: tuple(1 if i == j else 0 for i in range(NV)) + (0,)
    Em = []
    for j in range(NV):
        Em.append({tuple(k if i == j else 0 for i in range(NV)) + (0,): Q2(Fr(1, math.factorial(k))) for k in range(1, DEG + 1)})
    u = {}
    Fp = {}
    for j in range(NV):
        u = padd(u, pmul(tpoly(*T[j]), Em[j]))
        Fp = padd(Fp, pscale(Em[j], T[j][1]))
    invF = const(1); pw = const(1)
    for k in range(1, DEG + 1):
        pw = pmul(pw, pscale(u, Q2(-1)))
        invF = padd(invF, pw)
    g = pmul(Fp, invF)
    yD = tpoly(q0, L2); yy = padd(tpoly(V02, q0 * Q2(2)), {(0,) * NV + (2,): L2})
    gg = pmul(g, g)
    inner = padd(pscale(pmul(yD, g), Q2(-2) / L2), pscale(pmul(yy, gg), ONE / L2))
    sq = const(1); pw = const(1)
    for k in range(1, DEG + 1):
        pw = pmul(pw, inner)
        c = Fr(1)
        for i in range(k):
            c *= (Fr(1, 2) - i)
        c /= math.factorial(k)
        sq = padd(sq, pscale(pw, Q2(c)))
    arc = pmul(sq, invF)
    c5 = {k: v for k, v in arc.items() if sum(k[:NV]) == DEG}
    # L/u = 3 - 2 sqrt2 (L^2 = (3-2sqrt2)^2 (2+sqrt2))
    Lu = Q2(3, -2)
    assert Lu * Lu * Q2(2, 1) == L2
    c5 = {k: v * Lu for k, v in c5.items()}
    # substitution sigma_j = alpha_j X + beta_j Y + ell_j, ell in {a: log t, b: log(1-t), c: log(1+k0 t), d: log(1+k1 t)}
    sub = []
    for j in range(NV):
        kind = info[j][0]
        cst = info[j][1]
        al, be = unitvec(cst)
        if kind == 'a':
            slot = 0
        elif kind == 'b':
            slot = 1
        else:
            k = info[j][2]
            slot = 2 if k == kap0 else (3 if k == kap1 else None)
            assert slot is not None, k
        sub.append((al, be, slot))
    print('sigma_j = alpha X + beta Y + ell_slot:', sub)
    rows = {}
    for key, v in c5.items():
        exps, n = key[:NV], key[NV]
        # expand prod_j (al_j X + be_j Y + ell_j)^{e_j}
        terms = {(0, 0, 0, 0, 0, 0): Q2(1)}
        for j in range(NV):
            al, be, slot = sub[j]
            for _ in range(exps[j]):
                new = {}
                for mon, c in terms.items():
                    for add_idx, f in ((4, al), (5, be), (slot, Fr(1))):
                        if f == 0:
                            continue
                        m2 = list(mon); m2[add_idx] += 1; m2 = tuple(m2)
                        new[m2] = new.get(m2, Q2(0)) + c * Q2(f)
                terms = new
        for mon, c in terms.items():
            if c.iszero():
                continue
            rows.setdefault(mon, {})
            rows[mon][n] = rows[mon].get(n, Q2(0)) + c * v
    rows = {m: {n: c for n, c in d.items() if not c.iszero()} for m, d in rows.items()}
    rows = {m: d for m, d in rows.items() if d}
    pickle.dump(rows, open(os.path.join(os.environ['G9WORK'], 'own_weight5_rows.pkl'), 'wb'))
    parsed = pickle.load(open('<project>/research/08_theory/OCT5D_octagon_order5_depth/parsed_mons.pkl', 'rb'))
    theirs = {}
    for mon, cu in parsed:
        if sum(mon) != 5:
            continue
        theirs[tuple(mon)] = {n: Q2(Fr(str(p)), Fr(str(q))) for n, (p, q) in cu.items()}
    nbad = 0; ncmp = 0
    for k in set(rows) | set(theirs):
        A = rows.get(k, {}); B = theirs.get(k, {})
        for n in set(A) | set(B):
            ncmp += 1
            if not (A.get(n, Q2(0)) - B.get(n, Q2(0))).iszero():
                nbad += 1
                if nbad <= 10:
                    print('  mismatch', k, n, A.get(n), B.get(n))
    print('weight-5 rows: own %d, OCT5D table %d; %d coefficients compared, mismatches: %d' % (len(rows), len(theirs), ncmp, nbad))


if __name__ == '__main__':
    main()
