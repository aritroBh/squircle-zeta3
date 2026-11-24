fix = r'''
def pf_from_sympy(Rexpr):
    """sympy rational in v over Q(sqrt2) -> PF. Only simple real poles in {0,+-a,+-b}; (v^2+1)^m."""
    from math import comb
    Rexpr = sp.together(Rexpr)
    num, den = sp.fraction(Rexpr)
    num = sp.expand(num); den = sp.expand(den)
    npoly = sp.Poly(num, v); dpoly = sp.Poly(den, v)
    assert npoly.degree() < dpoly.degree(), "improper fraction"
    N = list(reversed([qsi_from_sympy_real(c) for c in npoly.all_coeffs()]))
    D = list(reversed([qsi_from_sympy_real(c) for c in dpoly.all_coeffs()]))
    def peval(poly, x):
        tot = ZERO
        for c in reversed(poly):
            tot = tot * x + c
        return tot
    def pdeflate(poly, x):
        n = len(poly) - 1
        q = [ZERO]*n
        carry = poly[n]
        q[n-1] = carry
        for t in range(n-1, 0, -1):
            carry = poly[t] + carry * x
            q[t-1] = carry
        rem = poly[0] + carry * x
        assert rem.is_zero(), ("deflation remainder", x, rem)
        return q
    def pmul(P1, P2):
        out = [ZERO]*(len(P1)+len(P2)-1)
        for a, ca in enumerate(P1):
            for b, cb in enumerate(P2):
                out[a+b] = out[a+b] + ca*cb
        return out
    P = pf_zero()
    Drem = list(D)
    roots_found = []
    for r in (ZERO, AQ, -AQ, BQ, -BQ):
        while len(Drem) >= 2 and peval(Drem, r).is_zero():
            Drem2 = pdeflate(Drem, r)
            if peval(Drem2, r).is_zero():
                raise ValueError(f"higher real pole at {r}")
            res = peval(N, r) / peval(Drem2, r)
            P[0][r] = P[0].get(r, ZERO) + res
            Drem = Drem2
            roots_found.append(r)
    for r, c in list(P[0].items()):
        if c.is_zero(): del P[0][r]
    m = (len(Drem) - 1)//2
    assert 2*m + 1 == len(Drem), ("odd leftover denominator", Drem)
    lead = Drem[-1]
    pw = [ONE]
    for _ in range(m):
        new = [ZERO]*(len(pw)+2)
        for t, c in enumerate(pw):
            new[t] = new[t] + c
            new[t+2] = new[t+2] + c
        pw = new
    for a1, b1 in zip(pw, Drem):
        assert (a1*lead - b1).is_zero(), ("denominator not lead*(v^2+1)^m", Drem)
    # full denominator D(v) = lead (v^2+1)^m prod_t (v - r_t)
    # at v = i+w: D = w^m * H(w),  H(w) = lead (w+2i)^m prod_t (w + i - r_t)
    H = [lead]
    for _ in range(m):
        H = pmul(H, [2*IQ, ONE])         # (w + 2i)
    for r in roots_found:
        H = pmul(H, [IQ - r, ONE])       # (w + i - r)
    # series of N(i+w) mod w^m
    Nw = []
    for s in range(m):
        tot = ZERO
        for t, c in enumerate(N):
            if t >= s:
                tot = tot + c * Fr(comb(t, s)) * (IQ ** (t-s))
        Nw.append(tot)
    # q = Nw / H mod w^m
    q = []
    for s in range(m):
        acc = Nw[s]
        for t in range(1, s+1):
            acc = acc - (H[t] if t < len(H) else ZERO) * q[s-t]
        q.append(acc / H[0])
    ipos = [q[m-j] for j in range(1, m+1)]
    while ipos and ipos[-1].is_zero(): ipos.pop()
    P[1] = ipos
    P[2] = [e.conj_i() for e in ipos]
    return P
'''
import re
src = open('fastred.py').read()
start = src.index('def pf_from_sympy')
end = src.index('def pf_to_sympy')
src = src[:start] + fix.strip() + '\n\n' + src[end:]
open('fastred.py','w').write(src)
print("patched")
