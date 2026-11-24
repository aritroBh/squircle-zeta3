"""amzv.py -- exact reduction of alternating MZVs of weight <= 4 (level 2), as iterated integrals
   Z(w) := G(w; 1),  w a word over {0, 1, -1}, first letter != 1, last letter != 0  (convergent),
G(a_1..a_n; z) = int_0^z dt/(t-a_1) G(a_2..a_n; t)  (a_1 outermost).
Relations used (all theorems):
  * shuffle products of convergent words (iterated integrals);
  * stuffle (quasi-shuffle) products of convergent nested sums, via
      G(0^{n1-1},a1,...,0^{nk-1},ak; 1) = (-1)^k Li_{n1..nk}(1/a1, a1/a2, ..., a_{k-1}/ak),
      Li_{s}(x) = sum_{m1>...>mk>=1} prod x_i^{m_i}/m_i^{s_i};
  * distribution (level 2):  sum_{eps in {+-1}^k} Li_s(eps_1 y_1,...) = 2^{k-|s|} Li_s(y_1^2,...)  (s_1 >= 2);
  * shuffle-regularisation  G(1;1) := 0  and the stuffle-regularisation  Li_1(1) := T  (formal), with the
    standard comparison at weights <= 4 not needed: we only use products in which every factor converges.
Basis constants: lam = log2, z2, z3, z4, li4h = Li4(1/2).  Products of lower-weight constants are kept as
monomials.  Identification of the free weight-4 element: G(-1,1,1,1;1) = -Li4(1/2) (duality t -> 1-t).
Everything exact (Fractions); a numerical self-check (mpmath, 40 digits) is printed by main().
"""
from fractions import Fraction as Fr
from collections import defaultdict
from functools import lru_cache
import itertools

CONSTS = ('lam', 'z2', 'z3', 'z4', 'li4h')
WT = {'lam': 1, 'z2': 2, 'z3': 3, 'z4': 4, 'li4h': 4}


def mono_mul(m1, m2):
    return tuple(a + b for a, b in zip(m1, m2))


ZERO_M = (0,) * len(CONSTS)


def unit(name):
    m = [0] * len(CONSTS)
    m[CONSTS.index(name)] = 1
    return tuple(m)


def padd(*ps, coeffs=None):
    out = defaultdict(Fr)
    for i, p in enumerate(ps):
        c = Fr(1) if coeffs is None else Fr(coeffs[i])
        for m, v in p.items():
            out[m] += c * v
    return {m: v for m, v in out.items() if v}


def pmul(p, q):
    out = defaultdict(Fr)
    for m1, a in p.items():
        for m2, b in q.items():
            out[mono_mul(m1, m2)] += a * b
    return {m: v for m, v in out.items() if v}


def pconst(c):
    return {ZERO_M: Fr(c)} if c else {}


@lru_cache(maxsize=None)
def shuffle_words(u, v):
    if not u:
        return ((v, 1),)
    if not v:
        return ((u, 1),)
    out = defaultdict(int)
    for w, c in shuffle_words(u[1:], v):
        out[(u[0],) + w] += c
    for w, c in shuffle_words(u, v[1:]):
        out[(v[0],) + w] += c
    return tuple(out.items())


def convergent(w):
    return len(w) > 0 and w[0] != 1 and w[-1] != 0


def words(n):
    return [w for w in itertools.product((0, 1, -1), repeat=n) if convergent(w)]


def to_sum(w):
    """G-word -> (sign, s tuple, x tuple) with G(w;1) = sign * Li_s(x)"""
    s, a = [], []
    cnt = 0
    for l in w:
        cnt += 1
        if l != 0:
            s.append(cnt)
            a.append(l)
            cnt = 0
    assert cnt == 0
    x = [Fr(1, a[0])] + [Fr(a[i - 1], a[i]) for i in range(1, len(a))]
    return (-1) ** len(a), tuple(s), tuple(int(v) for v in x)


def to_word(s, x):
    """Li_s(x) (x_i in {+-1}) -> (sign, G-word)"""
    a = []
    prod = 1
    for xi in x:
        prod *= xi
        a.append(prod)            # a1 = 1/x1 = x1 (x1 = +-1), a_i = a_{i-1}/x_i = a_{i-1}*x_i
    w = []
    for si, ai in zip(s, a):
        w += [0] * (si - 1) + [ai]
    return (-1) ** len(s), tuple(w)


def stuffle(A, B):
    """quasi-shuffle of (s,x) sequences: list of ((s,x), coeff)"""
    (sa, xa), (sb, xb) = A, B
    if not sa:
        return [((sb, xb), 1)]
    if not sb:
        return [((sa, xa), 1)]
    out = []
    for (s, x), c in stuffle((sa[1:], xa[1:]), (sb, xb)):
        out.append((((sa[0],) + s, (xa[0],) + x), c))
    for (s, x), c in stuffle((sa, xa), (sb[1:], xb[1:])):
        out.append((((sb[0],) + s, (xb[0],) + x), c))
    for (s, x), c in stuffle((sa[1:], xa[1:]), (sb[1:], xb[1:])):
        out.append((((sa[0] + sb[0],) + s, (xa[0] * xb[0],) + x), c))
    return out


class AMZV:
    def __init__(self, maxw=4):
        self.val = {}          # convergent word -> poly in constants
        self.maxw = maxw
        self._build()

    def value(self, w):
        w = tuple(w)
        if not w:
            return pconst(1)
        return self.val[w]

    def _build(self):
        # weight 1: G(-1;1) = log 2
        self.val[(-1,)] = {unit('lam'): Fr(1)}
        for n in range(2, self.maxw + 1):
            self._solve_weight(n)

    def _lower(self, w):
        return self.value(w)

    def _solve_weight(self, n):
        W = words(n)
        idx = {w: i for i, w in enumerate(W)}
        rels = []                       # list of (dict word->coef, poly rhs) meaning sum coef*Z(w) = rhs
        # shuffle relations: u, v convergent, |u|+|v| = n
        for k in range(1, n):
            for u in words(k):
                for v in words(n - k):
                    if u > v:
                        continue
                    lhs = defaultdict(Fr)
                    for w, c in shuffle_words(u, v):
                        lhs[w] += c
                    rels.append((dict(lhs), pmul(self.value(u), self.value(v))))
        # stuffle relations
        for k in range(1, n):
            for u in words(k):
                for v in words(n - k):
                    if u > v:
                        continue
                    su, s1, x1 = to_sum(u)
                    sv, s2, x2 = to_sum(v)
                    lhs = defaultdict(Fr)
                    for (s, x), c in stuffle((s1, x1), (s2, x2)):
                        sg, w = to_word(s, x)
                        lhs[w] += Fr(c * sg)
                    # Li(s1,x1) Li(s2,x2) = sum  ->  (su*Z(u))(sv*Z(v)) = sum sg*Z(w) ... careful with signs:
                    # Z(u) = su * Li(s1,x1)  =>  Li(s1,x1) = su*Z(u);  Li(s,x) = sg*Z(w)
                    rels.append((dict(lhs), pmul(pconst(su * sv), pmul(self.value(u), self.value(v)))))
        # distribution: sum_eps Li_s(eps) = 2^{k-|s|} Li_s(1,..,1), s_1 >= 2
        for k in range(1, n + 1):
            for s in itertools.product(range(1, n + 1), repeat=k):
                if sum(s) != n or s[0] < 2:
                    continue
                lhs = defaultdict(Fr)
                for eps in itertools.product((1, -1), repeat=k):
                    sg, w = to_word(s, eps)
                    lhs[w] += Fr(sg)
                sg1, w1 = to_word(s, (1,) * k)
                lhs[w1] -= Fr(2) ** (k - n) * sg1
                rels.append((dict(lhs), {}))
        # normalisation of free elements
        norm = []
        if n == 2:
            norm.append(({(0, 1): Fr(1)}, {unit('z2'): Fr(-1)}))
        if n == 3:
            norm.append(({(0, 0, 1): Fr(1)}, {unit('z3'): Fr(-1)}))
        if n == 4:
            norm.append(({(0, 0, 0, 1): Fr(1)}, {unit('z4'): Fr(-1)}))
            norm.append(({(-1, 1, 1, 1): Fr(1)}, {unit('li4h'): Fr(-1)}))
        rels += norm
        # exact elimination: unknowns = words of weight n; rhs in polynomial ring (vectorised over monomials)
        monos = sorted({m for _, r in rels for m in r} | {ZERO_M})
        nW = len(W)
        rows = []
        for lhs, rhs in rels:
            row = [Fr(0)] * nW
            for w, c in lhs.items():
                row[idx[w]] += c
            rows.append((row, dict(rhs)))
        # Gaussian elimination
        pivots = {}
        basis_rows = []
        for row, rhs in rows:
            row = row[:]
            rhs = dict(rhs)
            for (pc, prow, prhs) in basis_rows:
                if row[pc] != 0:
                    f = row[pc]
                    row = [a - f * b for a, b in zip(row, prow)]
                    rhs = padd(rhs, prhs, coeffs=[1, -f])
            nz = [j for j in range(nW) if row[j] != 0]
            if not nz:
                if rhs:
                    raise ValueError('inconsistent relation at weight %d: %s' % (n, rhs))
                continue
            pc = nz[0]
            f = row[pc]
            row = [a / f for a in row]
            rhs = {m: v / f for m, v in rhs.items()}
            # back-substitute into existing rows
            new_basis = []
            for (qc, qrow, qrhs) in basis_rows:
                if qrow[pc] != 0:
                    g = qrow[pc]
                    qrow = [a - g * b for a, b in zip(qrow, row)]
                    qrhs = padd(qrhs, rhs, coeffs=[1, -g])
                new_basis.append((qc, qrow, qrhs))
            basis_rows = new_basis + [(pc, row, rhs)]
        solved = {}
        free = set(range(nW)) - {pc for pc, _, _ in basis_rows}
        if free:
            raise ValueError('weight %d: %d words not determined, e.g. %s' % (n, len(free), [W[j] for j in sorted(free)][:5]))
        for pc, row, rhs in basis_rows:
            solved[W[pc]] = rhs
        self.val.update(solved)


def numeric_G(w, dps=40):
    import mpmath as mp
    mp.mp.dps = dps + 10
    sg, s, x = to_sum(tuple(w))
    # nested sum via mpmath polylog for depth 1, else direct series (alternating -> use nsum acceleration)
    if len(s) == 1:
        return sg * mp.polylog(s[0], x[0])
    # depth >= 2: iterated integral by recursion of HPL-type series at z -> 1 is slow; use mpmath quad on
    # the iterated integral representation recursively (adequate for 40 digits at weight <= 4)
    def G(word, z):
        if not word:
            return mp.mpf(1)
        a = word[0]
        if len(word) == 1:
            return mp.log(1 - z / a) if a != 0 else mp.log(z)
        return mp.quad(lambda t: G(word[1:], t) / (t - a), [0, z])
    return G(tuple(w), mp.mpf(1))


def poly_numeric(p, dps=40):
    import mpmath as mp
    mp.mp.dps = dps + 10
    vals = {'lam': mp.log(2), 'z2': mp.zeta(2), 'z3': mp.zeta(3), 'z4': mp.zeta(4), 'li4h': mp.polylog(4, 0.5)}
    tot = mp.mpf(0)
    for m, v in p.items():
        t = mp.mpf(v.numerator) / v.denominator
        for name, e in zip(CONSTS, m):
            t *= vals[name] ** e
        tot += t
    return tot


if __name__ == '__main__':
    import mpmath as mp
    A = AMZV(4)
    print('solved words:', len(A.val))
    worst = 0
    for w in sorted(A.val, key=lambda w: (len(w), w)):
        if len(w) <= 3 or w in [(-1, 0, 0, -1), (-1, 0, 0, 1), (0, 1, 0, -1), (-1, -1, -1, -1), (0, -1, 0, 1), (-1, 1, 0, 1)]:
            num = numeric_G(w, 30)
            ex = poly_numeric(A.val[w], 30)
            d = abs(num - ex)
            worst = max(worst, d)
    print('numerical check (sampled words, 30 digits): worst |G - exact| = %s' % mp.nstr(worst, 5))
