code = r'''
"""Fast exact IBP reducer for M_k integrands.

All rational functions that occur have:
  - simple real poles in {0, a=1+sqrt2, -a, b=sqrt2-1, -b}
  - poles of any order at v = +i and v = -i (conjugate principal parts)
  - zero polynomial part
We represent such a function in exact partial fractions with coefficients in
Q(sqrt2, i), stored as 4 Fractions (p0 + p1*sqrt2 + p2*i + p3*sqrt2*i).

Mirrors ibp_reduce.reduce_k exactly, but with no sympy apart/integrate.
"""
from fractions import Fraction as Fr
import sympy as sp

S2 = sp.sqrt(2)

# ---------- Q(sqrt2, i) arithmetic ----------
class QSI:
    __slots__ = ("p0", "p1", "p2", "p3")
    def __init__(s, p0=0, p1=0, p2=0, p3=0):
        s.p0 = Fr(p0); s.p1 = Fr(p1); s.p2 = Fr(p2); s.p3 = Fr(p3)
    @staticmethod
    def c(x):
        return x if isinstance(x, QSI) else QSI(x)
    def __add__(s, o):
        o = QSI.c(o); return QSI(s.p0+o.p0, s.p1+o.p1, s.p2+o.p2, s.p3+o.p3)
    __radd__ = __add__
    def __neg__(s):
        return QSI(-s.p0, -s.p1, -s.p2, -s.p3)
    def __sub__(s, o):
        return s + (-QSI.c(o))
    def __rsub__(s, o):
        return QSI.c(o) - s
    def scale(s, q):
        q = Fr(q); return QSI(s.p0*q, s.p1*q, s.p2*q, s.p3*q)
    def __mul__(s, o):
        o = QSI.c(o)
        a,b,c,d = s.p0,s.p1,s.p2,s.p3
        e,f,g,h = o.p0,o.p1,o.p2,o.p3
        return QSI(a*e + 2*b*f - c*g - 2*d*h,
                   a*f + b*e - c*h - d*g,
                   a*g + c*e + 2*b*h + 2*d*f,
                   a*h + d*e + b*g + c*f)
    __rmul__ = __mul__
    def conj_i(s):   # i -> -i
        return QSI(s.p0, s.p1, -s.p2, -s.p3)
    def conj_s(s):   # sqrt2 -> -sqrt2
        return QSI(s.p0, -s.p1, s.p2, -s.p3)
    def inv(s):
        # 1/s = conj_s / N, N in Q(i); then 1/N = conj_i/|N|^2
        cs = s.conj_s()
        N = s * cs                     # in Q(i): p1=p3=0
        assert N.p1 == 0 and N.p3 == 0
        n2 = N.p0*N.p0 + N.p2*N.p2     # |N|^2 in Q
        assert n2 != 0
        Ni = QSI(N.p0/n2, 0, -N.p2/n2, 0)
        return cs * Ni
    def __truediv__(s, o):
        return s * QSI.c(o).inv()
    def __rtruediv__(s, o):
        return QSI.c(o) / s
    def __pow__(s, k):
        r = QSI(1); base = s if k >= 0 else s.inv()
        for _ in range(abs(k)):
            r = r * base
        return r
    def is_zero(s):
        return s.p0 == 0 and s.p1 == 0 and s.p2 == 0 and s.p3 == 0
    def is_real(s):
        return s.p2 == 0 and s.p3 == 0
    def __eq__(s, o):
        o = QSI.c(o)
        return s.p0==o.p0 and s.p1==o.p1 and s.p2==o.p2 and s.p3==o.p3
    def __hash__(s):
        return hash((s.p0,s.p1,s.p2,s.p3))
    def __repr__(s):
        return f"QSI({s.p0},{s.p1},{s.p2},{s.p3})"
    def sympy(s):
        return sp.Rational(s.p0.numerator, s.p0.denominator) + sp.Rational(s.p1.numerator, s.p1.denominator)*sp.sqrt(2) \
             + sp.Rational(s.p2.numerator, s.p2.denominator)*sp.I + sp.Rational(s.p3.numerator, s.p3.denominator)*sp.sqrt(2)*sp.I

ZERO = QSI(0); ONE = QSI(1); IQ = QSI(0,0,1,0)
AQ = QSI(1,1,0,0)     # a = 1+sqrt2
BQ = QSI(-1,1,0,0)    # b = sqrt2-1

def qsi_from_sympy_real(expr):
    """exact sympy number in Q(sqrt2) -> QSI (real)"""
    e = sp.expand(sp.radsimp(expr))
    p0 = e.subs(S2, 0)
    p1 = sp.expand((e - p0)/S2)
    assert p0.free_symbols == set() and p1.free_symbols == set()
    p0 = sp.Rational(p0); p1 = sp.Rational(p1)
    return QSI(Fr(int(p0.p), int(p0.q)), Fr(int(p1.p), int(p1.q)))

# ---------- PF representation ----------
# PF = [realdict {QSI(real r): QSI(real c)}, ipos [e1..em], ineg [f1..fn], const QSI]
# meaning  sum c/(v-r) + sum_j ipos[j-1]/(v-i)^j + sum_j ineg[j-1]/(v+i)^j + const

def pf_zero():
    return [{}, [], [], ZERO]

def pf_copy(P):
    return [dict(P[0]), list(P[1]), list(P[2]), P[3]]

def pf_add(P, Q):
    R = pf_copy(P)
    for r, c in Q[0].items():
        R[0][r] = R[0].get(r, ZERO) + c
        if R[0][r].is_zero(): del R[0][r]
    def addlist(L1, L2):
        n = max(len(L1), len(L2))
        out = []
        for t in range(n):
            x = L1[t] if t < len(L1) else ZERO
            y = L2[t] if t < len(L2) else ZERO
            out.append(x + y)
        while out and out[-1].is_zero(): out.pop()
        return out
    R[1] = addlist(R[1], Q[1])
    R[2] = addlist(R[2], Q[2])
    R[3] = R[3] + Q[3]
    return R

def pf_neg(P):
    return [{r: -c for r, c in P[0].items()}, [-e for e in P[1]], [-f for f in P[2]], -P[3]]

def pf_sub(P, Q):
    return pf_add(P, pf_neg(Q))

def pf_scale(P, q):
    q = QSI.c(q)
    if q.is_zero(): return pf_zero()
    return [{r: c*q for r, c in P[0].items()}, [e*q for e in P[1]], [f*q for f in P[2]], P[3]*q]

def pf_check_reality(P, ctx=""):
    for r, c in P[0].items():
        assert r.is_real() and c.is_real(), (ctx, "real", r, c)
    for e, f in zip(P[1] + [ZERO]*max(0, len(P[2])-len(P[1])),
                    P[2] + [ZERO]*max(0, len(P[1])-len(P[2]))):
        assert e == f.conj_i(), (ctx, "conj", e, f)
    assert P[3].is_real(), (ctx, "const", P[3])

# fold: R(v) -> R(1/v)/v^2
def pf_fold(P):
    real, ipos, ineg, const = P
    assert const.is_zero(), "fold of nonzero constant"
    R = pf_zero()
    for r, c in real.items():
        if r.is_zero():
            R[0][r] = R[0].get(r, ZERO) + c
        else:
            # c/(v-r) -> c/v - c/(v - 1/r)
            R[0][ZERO] = R[0].get(ZERO, ZERO) + c
            rinv = r.inv()
            R[0][rinv] = R[0].get(rinv, ZERO) - c
    for r, c in list(R[0].items()):
        if c.is_zero(): del R[0][r]
    # +i poles -> -i poles:  e_j/(v-i)^j -> e_j v^{j-2}/((-i)^j (v+i)^j)
    def fold_side(coeffs, sign):
        # sign = +1: old pole at +i (v-i), maps to (v+i); sign = -1: old at -i maps to (v-i)
        # unit u with 1 - sign*i*v = u*(v + sign*i):  for sign=+1: 1-iv = -i(v+i), u=-i
        #                                             for sign=-1: 1+iv = i(v-i),  u= i
        u = QSI(0,0,-sign,0)   # -sign * i
        out_real0 = ZERO
        out = []   # new coeffs at (v + sign i)^{-order}
        def add_at(order, val):
            while len(out) < order: out.append(ZERO)
            out[order-1] = out[order-1] + val
        for j, e in enumerate(coeffs, start=1):
            if e.is_zero(): continue
            if j == 1:
                # e/[v*u*(v+sign i)]; 1/[v(v+sign i)] = (1/(sign i))(1/v - 1/(v+sign i))
                # total: e/u * (1/(sign i))(1/v - 1/(v+sign i))
                si = QSI(0,0,sign,0)
                c0 = e / u / si
                out_real0 = out_real0 + c0
                add_at(1, -c0)
            else:
                # e * v^{j-2}/(u^j (v+sign i)^j), v = w - sign i (w = v + sign i)
                si = QSI(0,0,sign,0)
                from math import comb
                for s in range(0, j-1):     # s = 0..j-2
                    coef = e * Fr(comb(j-2, s)) * ((-si) ** (j-2-s)) / (u ** j)
                    add_at(j - s, coef)
        return out_real0, out
    r0a, newneg = fold_side(ipos, +1)
    r0b, newpos = fold_side(ineg, -1)
    if not (r0a + r0b).is_zero():
        R[0][ZERO] = R[0].get(ZERO, ZERO) + r0a + r0b
    R[1] = newpos
    R[2] = newneg
    # trim
    for L in (R[1], R[2]):
        while L and L[-1].is_zero(): L.pop()
    for r, c in list(R[0].items()):
        if c.is_zero(): del R[0][r]
    return R

# split: simple part (real poles + order-1 at +-i) and higher part (order>=2 at +-i)
def pf_split(P):
    real, ipos, ineg, const = P
    simple = pf_zero(); higher = pf_zero()
    simple[0] = dict(real)
    if ipos: simple[1] = [ipos[0]]
    if ineg: simple[2] = [ineg[0]]
    higher[1] = [ZERO] + ipos[1:] if len(ipos) >= 2 else []
    higher[2] = [ZERO] + ineg[1:] if len(ineg) >= 2 else []
    for H in (higher[1], higher[2]):
        while H and H[-1].is_zero(): H.pop()
    return simple, higher, const

# integrate higher part: antiderivative rational (pure ipart)
def pf_integrate_higher(H):
    ipos, ineg = H[1], H[2]
    def integ(coeffs):
        out = []
        for j, e in enumerate(coeffs, start=1):
            if j == 1:
                assert e.is_zero(), "log term in higher antiderivative"
                out.append(ZERO)
            else:
                out.append(e.scale(Fr(-1, j-1)))
        while out and out[-1].is_zero(): out.pop()
        return out
    G = pf_zero()
    G[1] = integ(ipos)
    G[2] = integ(ineg)
    return G

# evaluate a pure-ipart function at v = b  (b real QSI)
def pf_eval_at(G, x):
    x = QSI.c(x)
    tot = ZERO
    for j, e in enumerate(G[1], start=1):
        tot = tot + e / ((x - IQ) ** j)
    for j, f in enumerate(G[2], start=1):
        tot = tot + f / ((x + IQ) ** j)
    tot = tot + G[3]
    for r, c in G[0].items():
        tot = tot + c / (x - r)
    return tot

# multiply (pure ipart + const) by DL = sum d_r/(v-r)
def pf_mul_DL(Gm, DL):
    realG, ipos, ineg, const = Gm
    R = pf_zero()
    for r, d in DL.items():
        if not const.is_zero():
            R[0][r] = R[0].get(r, ZERO) + const * d
        for coeffs, sgn in ((ipos, +1), (ineg, -1)):
            si = QSI(0,0,sgn,0)         # pole location  sgn*i
            den = r - si                # r - sgn*i
            for j, e in enumerate(coeffs, start=1):
                if e.is_zero(): continue
                ed = e * d
                R[0][r] = R[0].get(r, ZERO) + ed / (den ** j)
                for u in range(1, j+1):
                    val = -ed / (den ** (j - u + 1))
                    L = R[1] if sgn == +1 else R[2]
                    while len(L) < u: L.append(ZERO)
                    L[u-1] = L[u-1] + val
    for r, c in list(R[0].items()):
        if c.is_zero(): del R[0][r]
    for L in (R[1], R[2]):
        while L and L[-1].is_zero(): L.pop()
    return R

# ---------- conversion from sympy ----------
v = sp.Symbol('v', positive=True)

def pf_from_sympy(Rexpr):
    """sympy rational in v over Q(sqrt2) -> PF. Only simple real poles; (v^2+1)^m allowed."""
    Rexpr = sp.together(Rexpr)
    num, den = sp.fraction(Rexpr)
    num = sp.expand(num); den = sp.expand(den)
    npoly = sp.Poly(num, v); dpoly = sp.Poly(den, v)
    if npoly.degree() >= dpoly.degree():
        q, r = sp.div(npoly, dpoly)
        assert sp.expand(r.as_expr()) == 0 or True
        raise ValueError("improper fraction not expected")
    nc = [qsi_from_sympy_real(c) for c in npoly.all_coeffs()]  # high to low
    dc = [qsi_from_sympy_real(c) for c in dpoly.all_coeffs()]
    # polynomial coefficients as QSI lists low-to-high
    def hi2lo(L):
        return list(reversed(L))
    N = hi2lo(nc); D = hi2lo(dc)
    # find real roots among 0, +-a, +-b (simple) and multiplicity of +-i
    P = pf_zero()
    import copy
    Drem = list(D)
    def peval(poly, x):
        tot = ZERO
        for c in reversed(poly):
            tot = tot * x + c
        return tot
    def pdeflate(poly, x):
        # divide by (v - x), exact; returns quotient (low-to-high)
        n = len(poly) - 1
        q = [ZERO]*n
        carry = poly[n]
        q[n-1] = carry
        for t in range(n-1, 0, -1):
            carry = poly[t] + carry * x
            q[t-1] = carry
        # remainder check
        rem = poly[0] + carry * x
        assert rem.is_zero(), ("deflation remainder", x, rem)
        return q
    for r in (ZERO, AQ, -AQ, BQ, -BQ):
        while True:
            if len(Drem) >= 2 and peval(Drem, r).is_zero():
                # residue: N(r)/(Drem'(r) ...) -- simple pole assumed
                Drem2 = pdeflate(Drem, r)
                if peval(Drem2, r).is_zero():
                    raise ValueError(f"higher real pole at {r}")
                res = peval(N, r) / peval(Drem2, r)
                P[0][r] = P[0].get(r, ZERO) + res
                Drem = Drem2
            else:
                break
    # remaining denominator should be c*(v^2+1)^m
    m = (len(Drem) - 1)//2
    lead = Drem[-1]
    # check Drem = lead * (v^2+1)^m
    chk = [ONE, ZERO, ONE]
    pw = [ONE]
    for _ in range(m):
        pw = ([ZERO]+pw)
        pw2 = list(pw)
        for t in range(len(pw2)): 
            pass
        # multiply current power by (v^2+1): shift by 2 and add
        base = pw[:-0] if False else None
        break
    # simpler: build (v^2+1)^m by repeated multiplication
    pw = [ONE]
    for _ in range(m):
        new = [ZERO]*(len(pw)+2)
        for t, c in enumerate(pw):
            new[t] = new[t] + c
            new[t+2] = new[t+2] + c
        pw = new
    assert len(pw) == len(Drem)
    for a1, b1 in zip(pw, Drem):
        assert (a1*lead - b1).is_zero(), ("denominator not (v^2+1)^m", Drem)
    # series at v = i: N(i+w)/(lead*( (i+w)^2+1 )^m) = N(i+w)/(lead * w^m (w+2i)^m)
    # H(w) = lead*(w+2i)^m ; series q = N(i+w)/H mod w^m ; e_j = q_{m-j}
    Nw = []
    for s in range(m):   # N(i+w) coefficients up to w^{m-1}
        tot = ZERO
        for t, c in enumerate(N):
            if t >= s:
                from math import comb
                tot = tot + c * Fr(comb(t, s)) * (IQ ** (t-s))
        Nw.append(tot)
    # H(w) = lead*(w+2i)^m coefficients up to m-1
    from math import comb
    Hw = []
    for s in range(m):
        Hw.append(lead * Fr(comb(m, s)) * ((IQ*2) ** (m - s)))
    # series division q = Nw / Hw mod w^m
    q = []
    for s in range(m):
        acc = Nw[s] if s < len(Nw) else ZERO
        for t in range(1, s+1):
            acc = acc - Hw[t] * q[s-t]
        q.append(acc / Hw[0])
    ipos = []
    for j in range(1, m+1):
        ipos.append(q[m-j])
    while ipos and ipos[-1].is_zero(): ipos.pop()
    P[1] = ipos
    P[2] = [e.conj_i() for e in ipos]
    return P

def pf_to_sympy(P):
    real, ipos, ineg, const = P
    out = const.sympy()
    for r, c in real.items():
        out = out + c.sympy()/(v - r.sympy())
    for j, e in enumerate(ipos, start=1):
        out = out + e.sympy()/(v - sp.I)**j
    for j, f in enumerate(ineg, start=1):
        out = out + f.sympy()/(v + sp.I)**j
    return out

# ---------- the reducer ----------
DL1 = {BQ: ONE, -AQ: ONE, ZERO: -ONE}     # d log t / dv
DL2 = {AQ: ONE, -BQ: ONE, ZERO: -ONE}     # d log(1-t) / dv

def fold_key(key, R):
    i, j = key
    if i >= j:
        return key, R
    return (j, i), pf_fold(R)

def reduce_k_fast(k, integrand_pf=None):
    """integrand_pf: optional dict {(i,j): PF}. Returns dict {(i,j): PF} reduced."""
    if integrand_pf is None:
        from mk_integrands import v_integrand
        raw = v_integrand(k)
        integrand_pf = {key: pf_from_sympy(R) for key, R in raw.items()}
    work = {}
    for key, R in integrand_pf.items():
        kk, RR = fold_key(key, R)
        work[kk] = pf_add(work.get(kk, pf_zero()), RR)
    final = {}
    for deg in range(k, -1, -1):
        for key in sorted([kk for kk in work if sum(kk) == deg], reverse=True):
            R = work.pop(key)
            if deg == 0:
                final[key] = R
                continue
            simple, higher, const = pf_split(R)
            assert const.is_zero(), f"polynomial part at {key}"
            if simple[0] or simple[1] or simple[2]:
                final[key] = pf_add(final.get(key, pf_zero()), simple)
            if not higher[1] and not higher[2]:
                continue
            G = pf_integrate_higher(higher)
            # check: derivative matches higher by construction
            C = pf_eval_at(G, BQ)
            Gm = [G[0], G[1], G[2], -C]
            i, j = key
            if i >= 1:
                kk, RR = fold_key((i-1, j), pf_scale(pf_mul_DL(Gm, DL1), Fr(-i)))
                work[kk] = pf_add(work.get(kk, pf_zero()), RR)
            if j >= 1:
                kk, RR = fold_key((i, j-1), pf_scale(pf_mul_DL(Gm, DL2), Fr(-j)))
                work[kk] = pf_add(work.get(kk, pf_zero()), RR)
    for key in final:
        pf_check_reality(final[key], ctx=str(key))
    return final
'''
with open('fastred.py','w') as f:
    f.write(code)
print("written", len(code))
