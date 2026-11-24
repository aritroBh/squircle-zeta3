"""cwbuild.py -- build folded integrand PF dict for level k, and kernel->value utilities."""
from cwintegrand import CornerIntegrand, p_comp_1mt

def build_integrand_pf(ci, red, k):
    out = {}
    for i in range(k, (k-1)//2, -1):
        j = k - i
        Ti = ci.T_poly(k, i)
        num, den = ci.sqrt_piece_nd(k, Ti)
        P = red.pf_from_poly(num, den)
        if i != j:
            Tj = p_comp_1mt(ci.T_poly(k, j))
            num, den = ci.sqrt_piece_nd(k, Tj)
            P = red.pf_add(P, red.pf_from_poly(num, den))
        if j == 0:
            numA, denA = ci.AB_nd(k)
            PA = red.pf_from_poly(numA, denA)
            P = red.pf_add(P, PA)
            P = red.pf_add(P, PA)
        out[(i, j)] = P
    return out


# ---- plain-tuple PF storage (picklable) ----
def q_plain(q):
    if hasattr(q, 'a'):   # Quad
        return ((q.a.numerator, q.a.denominator), (q.b.numerator, q.b.denominator))
    return tuple((x.numerator, x.denominator) for x in (q.p0, q.p1, q.p2, q.p3))

def pf_plain(P):
    return {'real': {(q_plain(r)[0] if hasattr(r,'a') else q_plain(r)): q_plain(c) for r, c in P[0].items()},
            'ppos': [q_plain(e) for e in P[1]], 'pneg': [q_plain(e) for e in P[2]], 'const': q_plain(P[3])}
