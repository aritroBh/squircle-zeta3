"""Exact weight-<=3 bookkeeping over Q(sqrt2) and exact evaluation of
        int_{v0}^{v1} log|v-x| log|v-y| dv/(v-z)
for real S-unit letters x, y, z of Q(sqrt2) (and v0, v1 in Q(sqrt2)), in terms of Li3, Li2*log, log^3,
followed by reduction of all Li2/Li3 values with classical functional equations only.

Polynomial ring: keys are sorted tuples of atoms
    'lam' = log 2,  'L' = log(1+sqrt2),  'pi2' = pi^2,  'z3' = zeta(3),  ('Li2', z),  ('Li3', z)
with coefficients in Q(sqrt2) (class QS2).

Li3 identities used (real forms, each one is checked numerically when applied):
  (I3) inversion   Li3(z) - Li3(1/z) = -(pi^2/6) log(-z) - (1/6) log^3(-z),                     z < 0
  (T3) three-term  Li3(z) + Li3(1-z) + Li3(1-1/z) = z3 + log^3(z)/6 + (pi^2/6) log z - (1/2) log^2(z) log(1-z),  0<z<1
  (D3) duplication Li3(w^2) = 4 Li3(w) + 4 Li3(-w),                                              |w| <= 1
  (V3) Li3(0)=0, Li3(1)=z3, Li3(-1)=-(3/4) z3, Li3(1/2) = (7/8) z3 - (pi^2/12) log2 + (1/6) log^3 2.
Antiderivatives used (each one is verified here by exact symbolic differentiation):
  P1(s) = int log|s| ds/(s-1):  -Li2(1-s)  (s > 0);   log(-s) log(1-s) + Li2(s)  (s < 0)
  P2(s) = int log^2|s| ds/(s-1):
        0<s<1 :  log^2 s log(1-s) + 2 log s Li2(s) - 2 Li3(s)
        s<0   :  log^2(-s) log(1-s) + 2 log(-s) Li2(s) - 2 Li3(s)
        s>1   :  with r = 1/s:  -log^3 r/3 + log^2 r log(1-r) + 2 log r Li2(r) - 2 Li3(r)
"""
from fractions import Fraction as Fr
import mpmath as mp
import sympy as sp
from qsqrt2 import QS2, SQ2, A, B, sunit_log
import dilog_exact as DE

ONE = QS2(1)
ZERO = QS2(0)


def _akey(atom):
    order = {'lam': 0, 'L': 1, 'pi2': 2, 'z3': 3}
    if isinstance(atom, str):
        return (0, order[atom], 0.0)
    return (1 if atom[0] == 'Li2' else 2, 0, float(atom[1]))


class P:
    def __init__(self, d=None):
        self.d = {}
        for k, c in (d or {}).items():
            c = QS2.c(c)
            if c != ZERO:
                k = tuple(sorted(k, key=_akey))
                self.d[k] = self.d.get(k, ZERO) + c
        self.d = {k: c for k, c in self.d.items() if c != ZERO}

    @staticmethod
    def const(c):
        return P({(): c})

    def __add__(s, o):
        o = o if isinstance(o, P) else P.const(o)
        r = dict(s.d)
        for k, c in o.d.items():
            r[k] = r.get(k, ZERO) + c
        return P(r)
    __radd__ = __add__

    def __neg__(s):
        return P({k: -c for k, c in s.d.items()})

    def __sub__(s, o):
        return s + (-(o if isinstance(o, P) else P.const(o)))

    def __mul__(s, o):
        if not isinstance(o, P):
            o = QS2.c(o)
            return P({k: c * o for k, c in s.d.items()})
        r = {}
        for k1, c1 in s.d.items():
            for k2, c2 in o.d.items():
                k = tuple(sorted(k1 + k2, key=_akey))
                r[k] = r.get(k, ZERO) + c1 * c2
        return P(r)
    __rmul__ = __mul__

    def atoms(s, kind):
        out = set()
        for k in s.d:
            for a in k:
                if not isinstance(a, str) and a[0] == kind:
                    out.add(a[1])
        return out

    def subs_atom(s, atom, val):
        """replace atom by P val (polynomially)."""
        out = P()
        for k, c in s.d.items():
            term = P.const(c)
            for a in k:
                term = term * (val if a == atom else P({(a,): 1}))
            out = out + term
        return out

    def num(s):
        r2 = mp.sqrt(2)
        base = {'lam': mp.log(2), 'L': mp.log(1 + r2), 'pi2': mp.pi ** 2, 'z3': mp.zeta(3)}
        tot = mp.mpf(0)
        for k, c in s.d.items():
            x = mp.mpf(c.p.numerator) / c.p.denominator + r2 * mp.mpf(c.q.numerator) / c.q.denominator
            for a in k:
                if isinstance(a, str):
                    x *= base[a]
                else:
                    z = DE.float_exact(a[1])
                    val = mp.polylog(2 if a[0] == 'Li2' else 3, z)
                    assert abs(mp.im(val)) < mp.mpf(10) ** (-mp.mp.dps + 10), a
                    x *= mp.re(val)
            tot += x
        return tot

    def __repr__(s):
        if not s.d:
            return "0"
        parts = []
        for k in sorted(s.d, key=lambda k: (len(k), [_akey(a) for a in k])):
            nm = "*".join(a if isinstance(a, str) else f"{a[0]}({a[1]})" for a in k) or "1"
            parts.append(f"{s.d[k]}*{nm}")
        return " + ".join(parts)


LAM, LL, PI2, Z3 = P({('lam',): 1}), P({('L',): 1}), P({('pi2',): 1}), P({('z3',): 1})


def LOG(z):
    j, m = sunit_log(z)
    return LAM * Fr(j, 2) + LL * m


def LI2(z):
    z = QS2.c(z)
    if z == ZERO:
        return P()
    return P({(('Li2', z),): 1})


def LI3(z):
    z = QS2.c(z)
    if z == ZERO:
        return P()
    return P({(('Li3', z),): 1})


def from_W(w):
    """dilog_exact.W -> P"""
    out = P()
    for k, c in w.d.items():
        if k == ('1',):
            out = out + P.const(c)
        elif k[0] == 'Li2':
            out = out + LI2(k[1]) * c
        else:
            out = out + P({tuple(k): c})
    return out


# ------------------------------------------------------------------------------------------- antiderivatives
def P1(s):
    if s > ZERO:
        return -LI2(ONE - s)
    return LOG(-s) * LOG(ONE - s) + LI2(s)


def P2(s):
    if s == ZERO:
        return P()
    if s == ONE:                       # common limit from both sides
        return Z3 * (-2)
    if ZERO < s < ONE:
        ls = LOG(s)
        return ls * ls * LOG(ONE - s) + ls * LI2(s) * 2 - LI3(s) * 2
    if s < ZERO:
        ls = LOG(-s)
        return ls * ls * LOG(ONE - s) + ls * LI2(s) * 2 - LI3(s) * 2
    assert s > ONE
    r = ONE / s
    lr = LOG(r)
    return -(lr * lr * lr) * Fr(1, 3) + lr * lr * LOG(ONE - r) + lr * LI2(r) * 2 - LI3(r) * 2


def region(s):
    if s == ZERO:
        return 'zero'
    if s == ONE:
        return 'one'
    return 'neg' if s < ZERO else ('mid' if s < ONE else 'big')


def J1(c, w0, w1):
    """int_{w0}^{w1} log|w| dw/(w-c),  c != 0, path avoids 0 (may end there) and c."""
    s0, s1 = w0 / c, w1 / c
    _same_side(s0, s1)
    lc = LOG(c)
    return lc * (_logabs(s1 - ONE) - _logabs(s0 - ONE)) + _P1lim(s1, s0) - _P1lim(s0, s1)


def J2(c, w0, w1):
    """int_{w0}^{w1} log^2|w| dw/(w-c)."""
    s0, s1 = w0 / c, w1 / c
    _same_side(s0, s1)
    lc = LOG(c)
    return (lc * lc * (_logabs(s1 - ONE) - _logabs(s0 - ONE)) + lc * 2 * (_P1lim(s1, s0) - _P1lim(s0, s1))
            + P2(s1) - P2(s0))


def _logabs(x):
    # (log|s-1| at s = 1 only ever appears multiplied by log|c| = 0; return 0 there)
    return P() if x == ZERO else LOG(x)


def _same_side(s0, s1):
    r0, r1 = region(s0), region(s1)
    if 'one' in (r0, r1):
        other = r1 if r0 == 'one' else r0
        assert other in ('mid', 'big', 'one'), (s0, s1)
        return
    assert r0 == r1 or 'zero' in (r0, r1), (s0, s1)
    assert not ('zero' in (r0, r1) and 'big' in (r0, r1))


def _P1lim(s, other):
    """P1 at s; if s == 0 take the one-sided limit from the side of `other`."""
    if s == ZERO:
        return -PI2 * Fr(1, 6) if other > ZERO else P()
    return P1(s)


def J3pure(w0, w1):
    """int log^2|w| dw/w = log^3|w|/3 (w0, w1 != 0)."""
    l0, l1 = LOG(w0), LOG(w1)
    return (l1 * l1 * l1 - l0 * l0 * l0) * Fr(1, 3)


# -------------------------------------------------------------- the basic integrals on a subinterval
def int_logsq(x, z, v0, v1):
    """int_{v0}^{v1} log^2|v-x| dv/(v-z)."""
    if x == z:
        return J3pure(v0 - x, v1 - x)
    return J2(z - x, v0 - x, v1 - x)


def int_logratio_sq(x, y, z, v0, v1):
    """int_{v0}^{v1} log^2|(v-x)/(v-y)| dv/(v-z);  y must not be an endpoint."""
    assert y not in (v0, v1)
    w0, w1 = (v0 - x) / (v0 - y), (v1 - x) / (v1 - y)
    if z == y:            # d log|v-y| = -d log|w-1|
        return -J2(ONE, w0, w1)
    if z == x:            # d log|v-x| = d log|w| - d log|w-1|
        return J3pure(w0, w1) - J2(ONE, w0, w1)
    c1 = (x - z) / (y - z)
    return J2(c1, w0, w1) - J2(ONE, w0, w1)


def int_loglog(x, y, z, v0, v1):
    """int_{v0}^{v1} log|v-x| log|v-y| dv/(v-z)."""
    if x == y:
        return int_logsq(x, z, v0, v1)
    if y in (v0, v1):
        x, y = y, x
    return (int_logsq(x, z, v0, v1) + int_logsq(y, z, v0, v1) - int_logratio_sq(x, y, z, v0, v1)) * Fr(1, 2)


def int_log(x, z, v0, v1):
    """int_{v0}^{v1} log|v-x| dv/(v-z)."""
    if x == z:
        l0, l1 = LOG(v0 - x), LOG(v1 - x)
        return (l1 * l1 - l0 * l0) * Fr(1, 2)
    return J1(z - x, v0 - x, v1 - x)


def int_one(z, v0, v1):
    return LOG(v1 - z) - LOG(v0 - z)


def verify_antiderivatives():
    """exact symbolic check that P1', P2' are the claimed integrands on each region."""
    s, r = sp.symbols('s r', positive=True)
    res = []
    # region s>0 : P1 = -Li2(1-s)
    res.append(sp.simplify(sp.expand_func(sp.diff(-sp.polylog(2, 1 - s), s)) - sp.log(s) / (s - 1)))
    # region s<0 : s = -r
    P1n = sp.log(r) * sp.log(1 + r) + sp.polylog(2, -r)
    res.append(sp.simplify(sp.expand_func(sp.diff(P1n, r) * (-1)) - sp.log(r) / (-r - 1)))
    # 0<s<1
    P2m = sp.log(s) ** 2 * sp.log(1 - s) + 2 * sp.log(s) * sp.polylog(2, s) - 2 * sp.polylog(3, s)
    res.append(sp.simplify(sp.expand_func(sp.diff(P2m, s)) - sp.log(s) ** 2 / (s - 1)))
    # s<0
    P2n = sp.log(r) ** 2 * sp.log(1 + r) + 2 * sp.log(r) * sp.polylog(2, -r) - 2 * sp.polylog(3, -r)
    res.append(sp.simplify(sp.expand_func(sp.diff(P2n, r) * (-1)) - sp.log(r) ** 2 / (-r - 1)))
    # s>1, r = 1/s in (0,1): d/ds = -r^2 d/dr
    P2b = -sp.log(r) ** 3 / 3 + sp.log(r) ** 2 * sp.log(1 - r) + 2 * sp.log(r) * sp.polylog(2, r) - 2 * sp.polylog(3, r)
    res.append(sp.simplify(sp.expand_func(-r ** 2 * sp.diff(P2b, r)) - sp.log(r) ** 2 / (1 / r - 1)))
    return res


# ------------------------------------------------------------------------------------ reductions
CERT3 = []


def _check(name, lhs, rhs):
    with mp.workdps(60):
        err = lhs.num() - rhs.num()
    assert abs(err) < mp.mpf(10) ** -45, (name, lhs, rhs, err)


def li3_relations(z):
    """list of (name, P relation == 0) involving Li3(z), with all real arguments."""
    rels = []
    if z < ZERO:
        lz = LOG(-z)
        rels.append(("inversion (I3)", LI3(z) - LI3(ONE / z) + PI2 * lz * Fr(1, 6) + lz * lz * lz * Fr(1, 6)))
    if ZERO < z < ONE:
        for u in (z, ONE - z, ONE / (ONE - z)) if False else (z,):
            pass
    # three-term relations in which z occurs at any of the three slots
    cands = []
    if ZERO < z < ONE:
        cands += [z, ONE - z]
    if z < ZERO:
        cands.append(ONE / (ONE - z))          # 1 - 1/u = z  <=> u = 1/(1-z), in (0,1)
    for u in cands:
        if ZERO < u < ONE and u != ONE:
            try:
                lu, l1u = LOG(u), LOG(ONE - u)
            except ValueError:
                continue
            rel = (LI3(u) + LI3(ONE - u) + LI3(ONE - ONE / u) - Z3 - lu * lu * lu * Fr(1, 6)
                   - PI2 * lu * Fr(1, 6) + lu * lu * l1u * Fr(1, 2))
            rels.append((f"three-term (T3) at {u}", rel))
    # duplication
    if ZERO < z <= ONE:
        w = DE.sqrt_in_field(z) if z != ONE else None
        if w is not None and w != z:
            rels.append((f"duplication (D3) w={w}", LI3(z) - LI3(w) * 4 - LI3(-w) * 4))
    if -ONE <= z <= ONE and z != ZERO:
        rels.append((f"duplication (D3) w={z}", LI3(z * z) - LI3(z) * 4 - LI3(-z) * 4))
    return rels


LI3_VALUES = {
    ONE: Z3, -ONE: Z3 * Fr(-3, 4),
    QS2(Fr(1, 2)): Z3 * Fr(7, 8) - PI2 * LAM * Fr(1, 12) + LAM * LAM * LAM * Fr(1, 6),
}


def solve_li3(expr, basis, max_rounds=6, extra_relations=()):
    """Express all Li3 atoms of expr in terms of `basis` (list of arguments) + z3 + products,
    by Gaussian elimination over the closure of the classical relations.  Returns (P, used_relations)."""
    need = set(expr.atoms('Li3'))
    args = set(need)
    rels = {}
    frontier = set(need)
    for _ in range(max_rounds):
        new = set()
        for z in frontier:
            try:
                zrels = li3_relations(z)
            except ValueError:          # a non-S-unit log would be needed: skip this relation family
                zrels = []
            for nm, rel in zrels:
                key = (nm, tuple(sorted(float(a) for a in rel.atoms('Li3'))))
                if key in rels:
                    continue
                if any(abs(sunit_log(a)[0]) > 12 or abs(sunit_log(a)[1]) > 8 for a in rel.atoms('Li3')):
                    continue
                rels[key] = (nm, rel)
                for a in rel.atoms('Li3'):
                    if a not in args:
                        args.add(a); new.add(a)
        frontier = new
        if not frontier:
            break
    allrels = list(rels.values()) + list(extra_relations)
    for z, val in LI3_VALUES.items():
        allrels.append((f"value (V3) Li3({z})", LI3(z) - val))
    for nm, rel in allrels:
        _check(nm, rel, P())
    # Gaussian elimination: eliminate non-basis atoms
    basis = [QS2.c(b) for b in basis]
    order = [a for a in args if a not in basis]
    rows = [(nm, rel) for nm, rel in allrels]
    solved = {}
    used = []
    for a in sorted(order, key=lambda a: -abs(float(a)) if True else 0):
        atom = ('Li3', a)
        piv = None
        for idx, (nm, rel) in enumerate(rows):
            c = rel.d.get(((atom),), None)
            if c is not None and c != ZERO:
                piv = idx; break
        if piv is None:
            continue
        nm, rel = rows.pop(piv)
        c = rel.d[(atom,)]
        # atom = -(rel - c*atom)/c
        val = (rel - LI3(a) * c) * (QS2(-1) / c)
        solved[a] = val
        used.append(nm)
        rows = [(n2, r2.subs_atom(atom, val)) for n2, r2 in rows]
        for k in list(solved):
            solved[k] = solved[k].subs_atom(atom, val)
    out = expr
    for a in need:
        if a in basis:
            continue
        if a not in solved:
            raise RuntimeError(f"Li3({a}) not reducible with the relations at hand")
        out = out.subs_atom(('Li3', a), solved[a])
    # check that the result contains only basis Li3's
    leftover = [a for a in out.atoms('Li3') if a not in basis]
    if leftover:
        raise RuntimeError(f"leftover Li3 atoms {leftover}")
    return out, used


def reduce_li2_all(expr):
    """replace every Li2 atom by its reduction (dilog_exact)."""
    out = expr
    for z in list(expr.atoms('Li2')):
        val = from_W(DE.reduce_li2(z))
        out = out.subs_atom(('Li2', z), val)
    return out


# ------------------------------------------------------------------ proper linear-algebra reduction
def _li3_part(rel, atoms_index):
    vec = {}
    for k, c in rel.d.items():
        if len(k) == 1 and not isinstance(k[0], str) and k[0][0] == 'Li3':
            vec[atoms_index[k[0][1]]] = c
        else:
            for a in k:
                assert isinstance(a, str) or a[0] != 'Li3' or len(k) == 1, "Li3 inside a product?"
    return vec


def li3_closure(seed, max_rounds=6, extra_relations=()):
    args = set(seed)
    rels = {}
    frontier = set(seed)
    for _ in range(max_rounds):
        new = set()
        for z in frontier:
            try:
                zrels = li3_relations(z)
            except ValueError:
                zrels = []
            for nm, rel in zrels:
                key = (nm, tuple(sorted(float(a) for a in rel.atoms('Li3'))))
                if key in rels:
                    continue
                try:
                    if any(abs(sunit_log(a)[0]) > 12 or abs(sunit_log(a)[1]) > 8 for a in rel.atoms('Li3')):
                        continue
                except ValueError:
                    continue
                rels[key] = (nm, rel)
                for a in rel.atoms('Li3'):
                    if a not in args:
                        args.add(a); new.add(a)
        frontier = new
        if not frontier:
            break
    allrels = list(rels.values()) + list(extra_relations)
    for z, val in LI3_VALUES.items():
        allrels.append((f"value (V3) Li3({z})", LI3(z) - val))
    return allrels


def solve_li3_linear(expr, basis, extra_relations=(), max_rounds=6, verbose=False):
    """Find multipliers mu_r (in Q(sqrt2)) with  expr - sum mu_r rel_r  free of non-basis Li3 atoms.
    Returns (reduced P, list of (mu, name)) or raises with the obstruction."""
    basis = [QS2.c(b) for b in basis]
    seed = set(expr.atoms('Li3')) | set(basis)
    for _, r in extra_relations:
        seed |= r.atoms('Li3')
    allrels = li3_closure(seed, max_rounds, extra_relations)
    for nm, rel in allrels:
        _check(nm, rel, P())
    atoms = sorted({a for _, r in allrels for a in r.atoms('Li3')} | expr.atoms('Li3'), key=float)
    nonb = [a for a in atoms if a not in basis]
    idx = {a: i for i, a in enumerate(nonb)}
    # matrix columns = relations; rows = non-basis atoms; solve  M mu = t
    def vec_of(p):
        v = [ZERO] * len(nonb)
        for k, c in p.d.items():
            if len(k) == 1 and not isinstance(k[0], str) and k[0][0] == 'Li3' and k[0][1] in idx:
                v[idx[k[0][1]]] = v[idx[k[0][1]]] + c
        return v
    cols = [vec_of(r) for _, r in allrels]
    target = vec_of(expr)
    n, m = len(nonb), len(cols)
    # augmented matrix rows: n rows, m+1 columns
    Mx = [[cols[j][i] for j in range(m)] + [target[i]] for i in range(n)]
    piv_cols = []
    row = 0
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
        piv_cols.append(col)
        row += 1
        if row == n:
            break
    # consistency
    bad = [i for i in range(row, n) if Mx[i][m] != ZERO]
    if bad:
        raise RuntimeError("not in the span of the classical relations; obstruction rows: "
                           + "; ".join(str({nonb[j]: 0 for j in range(n)}) for _ in bad[:1])
                           + f"  (rank {row}, {n} non-basis atoms, inconsistent rows {len(bad)})")
    mu = [ZERO] * m
    for r, col in enumerate(piv_cols):
        mu[col] = Mx[r][m]
    red = expr
    used = []
    for j, (nm, rel) in enumerate(allrels):
        if mu[j] != ZERO:
            red = red - rel * mu[j]
            used.append((mu[j], nm))
    leftover = [a for a in red.atoms('Li3') if a not in basis]
    assert not leftover, leftover
    return red, used, (row, n)


# ------------------------------------------------- integrals with the pole dv/(v-e) at an endpoint e
def int_reg1(x, e, v0, v1):
    """int_{v0}^{v1} log|(v-x)/(e-x)| dv/(v-e),  e in {v0, v1}, x != e."""
    assert e in (v0, v1) and x != e
    w0, w1 = (v0 - x) / (e - x), (v1 - x) / (e - x)
    return J1(ONE, w0, w1)


def int_reg2(x, y, e, v0, v1):
    """int_{v0}^{v1} log|(v-x)/(e-x)| log|(v-y)/(e-y)| dv/(v-e),  e in {v0, v1}, x, y != e."""
    assert e in (v0, v1) and x != e and y != e
    def sq(xx):
        w0, w1 = (v0 - xx) / (e - xx), (v1 - xx) / (e - xx)
        return J2(ONE, w0, w1)
    if x == y:
        return sq(x)
    # w = (v-x)(e-y)/((v-y)(e-x)),  w(e) = 1,  w(inf) = (e-y)/(e-x);  d log|v-e| = dw/(w-1) - dw/(w-winf)
    assert y not in (v0, v1) or x not in (v0, v1)
    if y in (v0, v1):
        x, y = y, x
    k = (e - y) / (e - x)
    w0, w1 = (v0 - x) / (v0 - y) * k, (v1 - x) / (v1 - y) * k
    ratio = J2(ONE, w0, w1) - J2(k, w0, w1)
    return (sq(x) + sq(y) - ratio) * Fr(1, 2)
