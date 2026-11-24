"""Exact weight-2 bookkeeping over Q(sqrt2): logs of S-units, formal Li2 symbols, and a reducer that
rewrites every Li2 at a real S-unit argument into the basis
        {pi^2, log^2 2, log2*L, L^2, Li2(b)},   b = sqrt2 - 1,  L = log(1+sqrt2),
using ONLY the classical identities listed below.  Each rewrite is recorded as a human-readable
certificate line and is also checked numerically (the numerical check is a sanity check, not part
of the proof).

Identities (real-variable forms; Li2 real-valued on (-inf, 1]):
  (R)  reflection   Li2(z) + Li2(1-z) = pi^2/6 - log z log(1-z),             0 < z < 1
  (I)  inversion    Li2(z) + Li2(1/z) = -pi^2/6 - (1/2) log^2(-z),           z < 0
  (La) Landen       Li2(z) + Li2(z/(z-1)) = -(1/2) log^2(1-z),               z < 1
  (D)  duplication  Li2(w^2) = 2 Li2(w) + 2 Li2(-w),                         |w| <= 1
  (X)  Landen's chi_2 value   Li2(b) - Li2(-b) = pi^2/8 - L^2/2   (proved in the note: the
       fixed point of chi_2((1-x)/(1+x)) + chi_2(x) = pi^2/8 + (1/2) log x log((1+x)/(1-x)) )
  (V)  Li2(0)=0, Li2(1)=pi^2/6, Li2(-1)=-pi^2/12, Li2(1/2)=pi^2/12-(1/2)log^2 2.
  (A)  Abel's five-term relation, real form, valid on the convex set {x<1, y<1, x+y<1}:
       Li2(x/(1-y)) + Li2(y/(1-x)) - Li2(xy/((1-x)(1-y))) = Li2(x) + Li2(y) + log(1-x) log(1-y);
       used once, at (x, y) = (-sqrt2, 1-sqrt2), where it reads
       Li2(-1) + Li2(-b^2) - Li2(b^2) = Li2(-sqrt2) + Li2(-b) + (1/2) log2 * L.
"""
from fractions import Fraction as Fr
from collections import deque
import mpmath as mp
from qsqrt2 import QS2, SQ2, A, B, sunit_log

mp.mp.dps = 60
ONE = QS2(1)


class W:
    """Linear combination over Q(sqrt2) of weight-<=2 monomials.
    keys: ('pi2',), ('lam','lam'), ('L','lam')... sorted tuples of 'lam','L'; ('Li2', z) formal;
          ('lam',), ('L',), ('1',) for lower weight."""

    def __init__(self, d=None):
        self.d = {}
        for k, v in (d or {}).items():
            v = QS2.c(v)
            if v != QS2(0):
                self.d[k] = v

    def __add__(s, o):
        r = dict(s.d)
        for k, v in o.d.items():
            r[k] = r.get(k, QS2(0)) + v
        return W(r)

    def __neg__(s):
        return W({k: -v for k, v in s.d.items()})

    def __sub__(s, o):
        return s + (-o)

    def scale(s, c):
        c = QS2.c(c)
        return W({k: c * v for k, v in s.d.items()})

    def __rmul__(s, c):
        return s.scale(c)

    def li2_keys(s):
        return [k for k in s.d if k[0] == 'Li2']

    def num(s):
        r2 = mp.sqrt(2); L = mp.log(1 + r2); lam = mp.log(2)
        val = {'lam': lam, 'L': L}
        tot = mp.mpf(0)
        for k, c in s.d.items():
            cv = mp.mpf(c.p.numerator) / c.p.denominator + r2 * mp.mpf(c.q.numerator) / c.q.denominator
            if k == ('pi2',):
                x = mp.pi ** 2
            elif k == ('1',):
                x = mp.mpf(1)
            elif k[0] == 'Li2':
                z = k[1]
                x = mp.polylog(2, float_exact(z))
                assert mp.im(x) == 0 or abs(mp.im(x)) < mp.mpf(10) ** -50
                x = mp.re(x)
            else:
                x = mp.mpf(1)
                for f in k:
                    x *= val[f]
            tot += cv * x
        return tot

    def __repr__(s):
        parts = []
        order = {('pi2',): 0, ('lam', 'lam'): 1, ('L', 'lam'): 2, ('L', 'L'): 3}
        for k in sorted(s.d, key=lambda k: (order.get(k, 9), str(k))):
            name = {('pi2',): 'pi^2', ('lam', 'lam'): 'log^2(2)', ('L', 'lam'): 'log2*L',
                    ('L', 'L'): 'L^2', ('lam',): 'log2', ('L',): 'L', ('1',): '1'}.get(k)
            if name is None:
                name = f"Li2({k[1]})" if k[0] == 'Li2' else '*'.join(k)
            parts.append(f"{s.d[k]}*{name}")
        return " + ".join(parts) if parts else "0"


def float_exact(z):
    return mp.mpf(z.p.numerator) / z.p.denominator + mp.sqrt(2) * mp.mpf(z.q.numerator) / z.q.denominator


def LOG(z):
    """log|z| for an S-unit z, as a W of weight 1."""
    j, m = sunit_log(z)
    return W({('lam',): Fr(j, 2), ('L',): m})


def mul(u, v):
    """product of two weight-1 W's (only 'lam','L' keys)."""
    r = {}
    for k1, c1 in u.d.items():
        for k2, c2 in v.d.items():
            k = tuple(sorted(k1 + k2))
            r[k] = r.get(k, QS2(0)) + c1 * c2
    return W(r)


PI2 = W({('pi2',): 1})


def LI2(z):
    z = QS2.c(z)
    return W({('Li2', z): 1})


# ----------------------------------------------------------------------------------------------
# reducer
CERT = []          # certificate lines


def _record(name, z, rhs):
    lhs = LI2(z)
    with mp.workdps(60):
        err = lhs.num() - rhs.num()
    assert abs(err) < mp.mpf(10) ** -45, (name, z, err)
    CERT.append(f"  Li2({z}) = {rhs}      [{name}]")


def rule_R(z):          # 0<z<1
    if QS2(0) < z < ONE:
        rhs = PI2.scale(Fr(1, 6)) - mul(LOG(z), LOG(ONE - z)) - LI2(ONE - z)
        return 'reflection (R)', rhs
    return None


def rule_I(z):          # z<0
    if z < QS2(0):
        lz = LOG(-z)
        rhs = PI2.scale(Fr(-1, 6)) - mul(lz, lz).scale(Fr(1, 2)) - LI2(ONE / z)
        return 'inversion (I)', rhs
    return None


def rule_La(z):         # z<1, z != 0
    if z < ONE and z != QS2(0):
        l1 = LOG(ONE - z)
        rhs = -mul(l1, l1).scale(Fr(1, 2)) - LI2(z / (z - ONE))
        return 'Landen (La)', rhs
    return None


BASE = {
    QS2(0): W(),
    ONE: PI2.scale(Fr(1, 6)),
    QS2(-1): PI2.scale(Fr(-1, 12)),
    QS2(Fr(1, 2)): PI2.scale(Fr(1, 12)) - mul(LOG(QS2(2)), LOG(QS2(2))).scale(Fr(1, 2)),
    B: LI2(B),
    -B: LI2(B) - PI2.scale(Fr(1, 8)) + mul(LOG(A), LOG(A)).scale(Fr(1, 2)),
}
BASE_NAME = {QS2(0): 'value (V)', ONE: 'value (V)', QS2(-1): 'value (V)', QS2(Fr(1, 2)): 'value (V)',
             -B: "Landen's chi_2 value (X)"}


def sqrt_in_field(z):
    """w in Q(sqrt2) with w^2 = z and w>0, or None."""
    # (x + y sqrt2)^2 = x^2 + 2y^2 + 2xy sqrt2
    import sympy as sp
    x, y = sp.symbols('x y')
    sols = sp.solve([x ** 2 + 2 * y ** 2 - sp.Rational(z.p.numerator, z.p.denominator),
                     2 * x * y - sp.Rational(z.q.numerator, z.q.denominator)], [x, y], dict=True)
    for s in sols:
        if s[x].is_rational and s[y].is_rational:
            w = QS2(Fr(int(sp.numer(s[x])), int(sp.denom(s[x]))), Fr(int(sp.numer(s[y])), int(sp.denom(s[y]))))
            if w > QS2(0):
                return w
    return None


MEMO = {}
ABEL = {-B * B: (-SQ2, -B)}      # argument y/(1-x)  ->  (x, y)


def abel_rhs(z):
    x, y = ABEL[z]
    assert x < ONE and y < ONE and x + y < ONE          # inside the domain of validity
    assert y / (ONE - x) == z
    rhs = (LI2(x) + LI2(y) + mul(LOG(ONE - x), LOG(ONE - y)) - LI2(x / (ONE - y))
           + LI2(x * y / ((ONE - x) * (ONE - y))))
    return rhs


def reduce_li2(z, depth=0):
    """Return W in the basis, recording a certificate."""
    z = QS2.c(z)
    if z in MEMO:
        return MEMO[z]
    if z in BASE:
        if z in BASE_NAME and z != B:
            _record(BASE_NAME[z], z, BASE[z])
        MEMO[z] = BASE[z]
        return BASE[z]
    # BFS over single-term rules to reach a BASE argument or a field square
    prev = {z: None}
    dq = deque([z])
    target = None
    while dq:
        u = dq.popleft()
        if u in BASE or u in MEMO or u in ABEL or (u != z and (ONE >= u > QS2(0)) and sqrt_in_field(u) is not None
                                       and sqrt_in_field(u) != u):
            target = u
            break
        for rule in (rule_La, rule_I, rule_R):
            r = rule(u)
            if r is None:
                continue
            nxt = r[1].li2_keys()[0][1]
            if nxt not in prev and nxt <= ONE:
                prev[nxt] = (u, r)
                dq.append(nxt)
    if target is None:
        # z itself a square?
        w = sqrt_in_field(z) if QS2(0) < z <= ONE else None
        if w is None:
            raise RuntimeError(f"cannot reduce Li2({z})")
        target = z
    # resolve target value
    if target in MEMO:
        val = MEMO[target]
    elif target in BASE:
        val = reduce_li2(target)
    elif target in ABEL:
        rhs = abel_rhs(target)
        _record(f"Abel five-term (A) at (x,y)={ABEL[target]}", target, rhs)
        val = reduce_all(rhs)
        MEMO[target] = val
    else:
        w = sqrt_in_field(target)
        rhs = LI2(w).scale(2) + LI2(-w).scale(2)
        _record('duplication (D)', target, rhs)
        val = reduce_li2(w).scale(2) + reduce_li2(-w).scale(2)
        MEMO[target] = val
    # walk back the chain: Li2(u) = E - Li2(next)
    u = target
    while prev[u] is not None:
        par, (name, rhs) = prev[u]
        _record(name, par, rhs)
        elem = rhs + LI2(u)                   # rhs = E - Li2(u)  =>  E = rhs + Li2(u)
        val = elem - val
        MEMO[par] = val
        u = par
    return MEMO[z]


def reduce_all(w):
    """Substitute every formal Li2 in w by its reduction."""
    out = W({k: v for k, v in w.d.items() if k[0] != 'Li2'})
    for k in w.li2_keys():
        out = out + reduce_li2(k[1]).scale(w.d[k])
    return out
