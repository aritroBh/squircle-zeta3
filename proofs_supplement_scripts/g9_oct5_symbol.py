# Copied from research/08_theory/REF_HEX5_referee/own_symbol.py (code written by the referee of the notes, not by the authors of the notes).
# Changes: the default input is the table written by g9_oct5_rows.py in the directory $G9WORK.
"""Referee's own symbol computation for prop:octagon5depth (no project code).

Weight-five part of P^{(5)}(1)/u for the regular octagon = 960 * sum over rows (a,b,c,d,e,f) of
sum_n coef_n * log2^e * log(1+sqrt2)^f * int_0^1 t^n log^a t log^b(1-t) log^c(1-t/p) log^d(1-t/q) dt,
with p = 2+sqrt2 = -1/kappa0, q = -(1+sqrt2) = -1/kappa1 (rows: weight-5 rows of the integrand; the referee
checked them against his own leading-log expansion in oct_weight5_integrand.py).
Top weight: log(1-t/a) = G(a;t), products by shuffle, and int_0^1 t^n G(w;t) dt has top-weight part
(1 - w_1^{n+1})/(n+1) G(w;1) (integration by parts; w_1 = 1 gives 0).
Symbol of G(w;1) = I(0; w; 1) by Goncharov's formula
S(I(a_0; a_1..a_n; a_{n+1})) = sum_i S(I(..^a_i..)) (x) (a_i - a_{i+1}) - S(I(..^a_i..)) (x) (a_i - a_{i-1}),
zero entries dropped (shuffle regularisation).  Every entry is +-2^alpha (1+sqrt2)^beta; signs (torsion)
are dropped; the entry becomes alpha X + beta Y, X = log 2, Y = log(1+sqrt2).
Pairing (Section M4 of the paper): sum_w s(w) (<w1w2w3,A><w4w5,B> - <w3w4w5,A><w1w2,B>),
A in {[X,[X,Y]], [Y,[X,Y]]}, B = [X,Y].  Both orders of the symbol entries are reported.
"""
import pickle, itertools, sys
from fractions import Fraction as Fr
from functools import lru_cache


class Q2:
    """a + b sqrt2, a, b rational."""
    __slots__ = ('a', 'b')

    def __init__(self, a, b=0):
        self.a = Fr(a); self.b = Fr(b)

    def __add__(s, o):
        o = o if isinstance(o, Q2) else Q2(o); return Q2(s.a + o.a, s.b + o.b)
    __radd__ = __add__

    def __neg__(s):
        return Q2(-s.a, -s.b)

    def __sub__(s, o):
        return s + (-(o if isinstance(o, Q2) else Q2(o)))

    def __mul__(s, o):
        o = o if isinstance(o, Q2) else Q2(o)
        return Q2(s.a * o.a + 2 * s.b * o.b, s.a * o.b + s.b * o.a)
    __rmul__ = __mul__

    def __truediv__(s, o):
        o = o if isinstance(o, Q2) else Q2(o)
        n = o.a * o.a - 2 * o.b * o.b
        return s * Q2(o.a / n, -o.b / n)

    def __pow__(s, k):
        r = Q2(1)
        for _ in range(k):
            r = r * s
        return r

    def __eq__(s, o):
        o = o if isinstance(o, Q2) else Q2(o); return s.a == o.a and s.b == o.b

    def __hash__(s):
        return hash((s.a, s.b))

    def iszero(s):
        return s.a == 0 and s.b == 0

    def __repr__(s):
        return '(%s%s%s*sqrt2)' % (s.a, '+' if s.b > 0 else '-', abs(s.b)) if s.b else '%s' % s.a


ZERO, ONE = Q2(0), Q2(1)
P = Q2(2, 1)
QL = Q2(-1, -1)
LETTERS = {'0': ZERO, '1': ONE, 'p': P, 'q': QL}


def pow2half(a2):
    base = Q2(2) ** (abs(a2) // 2) * (Q2(0, 1) if abs(a2) % 2 else ONE)
    return base if a2 >= 0 else ONE / base


def unitvec(z):
    """z = +-2^alpha (1+sqrt2)^beta  ->  (alpha, beta)."""
    import math
    v = abs(float(z.a) + float(z.b) * math.sqrt(2))
    for a2 in range(-20, 21):
        for b in range(-20, 21):
            if abs(a2 / 2 * math.log(2) + b * math.log(1 + math.sqrt(2)) - math.log(v)) < 1e-11:
                u = Q2(1, 1) ** b if b >= 0 else ONE / Q2(1, 1) ** (-b)
                cand = pow2half(a2) * u
                if cand == z or cand == -z:
                    return (Fr(a2, 2), Fr(b))
    raise ValueError(z)


def sadd(d, k, c):
    v = d.get(k)
    v = c if v is None else v + c
    if v.iszero():
        d.pop(k, None)
    else:
        d[k] = v


@lru_cache(None)
def symI(seq):
    """symbol of I(a_0; a_1..a_n; a_{n+1}) as dict: tuple of entry-letters (pairs of Q2 diffs) -> int;
    entries are kept as unit vectors (alpha, beta)."""
    n = len(seq) - 2
    if n == 0:
        return {(): Q2(1)}
    out = {}
    for i in range(1, n + 1):
        sub = symI(seq[:i] + seq[i + 1:])
        for nb, sgn in ((seq[i + 1], 1), (seq[i - 1], -1)):
            dlt = seq[i] - nb
            if dlt.iszero():
                continue
            e = unitvec(dlt)
            if e == (0, 0):
                continue
            for w, c in sub.items():
                sadd(out, w + (e,), c * sgn)
    return out


def expand_words(sym):
    """entries (alpha,beta) -> words in letters 0 (=X) and 1 (=Y), multilinear."""
    out = {}
    for w, c in sym.items():
        choices = []
        for (al, be) in w:
            ch = []
            if al:
                ch.append((0, al))
            if be:
                ch.append((1, be))
            choices.append(ch)
        for combo in itertools.product(*choices):
            word = tuple(x for x, _ in combo)
            coef = c
            for _, f in combo:
                coef = coef * Q2(f)
            sadd(out, word, coef)
    return out


@lru_cache(None)
def symG_words(w):
    if not w:
        return {(): Q2(1)}
    """G(w_1..w_n;1) = int_0^1 dt/(t-w_1) G(w_2..;t) = I(0; w_n, ..., w_1; 1)."""
    seq = (ZERO,) + tuple(LETTERS[x] for x in reversed(w)) + (ONE,)
    return expand_words(symI(seq))


def shuffle_words(u, v):
    if not u:
        return {v: 1}
    if not v:
        return {u: 1}
    out = {}
    for w, c in shuffle_words(u[1:], v).items():
        k = (u[0],) + w; out[k] = out.get(k, 0) + c
    for w, c in shuffle_words(u, v[1:]).items():
        k = (v[0],) + w; out[k] = out.get(k, 0) + c
    return out


def shuffle_sym(s1, s2):
    out = {}
    for u, a in s1.items():
        for v, b in s2.items():
            for w, c in shuffle_words(u, v).items():
                sadd(out, w, a * b * Q2(c))
    return out


def letter_shuffles(a, b, c, d):
    """distinct words of 0^a sh 1^b sh p^c sh q^d with multiplicity (shuffle of the constant words)."""
    base = {(): 1}
    for let, k in (('0', a), ('1', b), ('p', c), ('q', d)):
        new = {}
        for u, cu in base.items():
            for w, cw in shuffle_words(u, (let,) * k).items():
                new[w] = new.get(w, 0) + cu * cw
        base = new
    return base


def lie():
    X = {(0,): Q2(1)}; Y = {(1,): Q2(1)}

    def concat(s, t):
        o = {}
        for u, a in s.items():
            for v, b in t.items():
                sadd(o, u + v, a * b)
        return o

    def br(s, t):
        o = dict(concat(s, t))
        for k, v in concat(t, s).items():
            sadd(o, k, -v)
        return o
    B = br(X, Y)
    return [br(X, B), br(Y, B)], B


def pairing(s, reverse=False):
    As, B = lie()
    if reverse:
        s = {w[::-1]: c for w, c in s.items()}
    res = []
    for A in As:
        tot = Q2(0)
        for w, c in s.items():
            t1 = A.get(w[:3], ZERO) * B.get(w[3:], ZERO)
            t2 = A.get(w[2:], ZERO) * B.get(w[:2], ZERO)
            tot = tot + c * (t1 - t2)
        res.append(tot)
    return res


def li_n_symbol(n, x):
    """Li_n(x) = -G(0^{n-1}, 1/x; 1): symbol -(1-x) (x) x (x) ... (x) x (first entry 1-x)."""
    e1 = unitvec(ONE - x); e2 = unitvec(x)
    return expand_words({(e1,) + (e2,) * (n - 1): Q2(-1)})


if __name__ == '__main__':
    import math
    # controls
    # (i) G(0,0,0,0,p;1) = -Li5(1/p)
    s = symG_words(('0', '0', '0', '0', 'p'))
    ref = li_n_symbol(5, ONE / P)
    diff = dict(s)
    for k, v in ref.items():
        sadd(diff, k, v)
    print('control G(0^4,p;1) + Li5(1/p) symbol = 0:', not diff)
    s2 = symG_words(('p', 'q')); s3 = symG_words(('q', 'p'))
    sh = shuffle_sym(symG_words(('p',)), symG_words(('q',)))
    d2 = dict(sh)
    for k, v in list(s2.items()) + list(s3.items()):
        sadd(d2, k, -v)
    print('control G(p)G(q) = G(p,q)+G(q,p) symbol:', not d2)
    b = Q2(-1, 1)
    for x in [b, -b, -(b * b), Q2(0, Fr(1, 2)), Q2(Fr(1, 2)), ONE - Q2(0, Fr(1, 2)), Q2(2, -1)]:
        pr = pairing(li_n_symbol(5, x)); prr = pairing(li_n_symbol(5, x), True)
        assert all(z.iszero() for z in pr + prr), (x, pr, prr)
    print('control: Li5 at 7 S-units pairs to 0 (both orders)')
    import random
    random.seed(3)
    okp = True
    for _ in range(50):
        k = random.choice([1, 2])
        w1 = tuple(random.choice('01pq') for _ in range(k)); w2 = tuple(random.choice('01pq') for _ in range(5 - k))
        if w1[0] == '1' or w2[0] == '1':
            continue
        pr = pairing(shuffle_sym(symG_words(w1), symG_words(w2)))
        okp &= all(z.iszero() for z in pr)
    print('control: shuffle products of G words pair to 0:', okp)
    # Li_{3,2}(x,1) = G(0,0,1/x,0,1/x;1)?  Li_{s1,s2}(x,y) = (-1)^2 G(0^{s1-1}, 1/x, 0^{s2-1}, 1/(xy); 1)
    # weight-5 rows of the octagon integrand
    import os
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.environ['G9WORK'], 'own_weight5_rows.pkl')
    rows = pickle.load(open(src, 'rb'))
    if isinstance(rows, dict):
        rows = list(rows.items())
    blocks = {'local': {}, 'single-p': {}, 'single-q': {}, 'mixed': {}}
    for mon, cu in rows:
        if sum(mon) != 5:
            continue
        a, b_, c, d, e, f = mon
        cls = 'mixed' if (c and d) else (('single-p' if c else 'single-q') if (c or d) else 'local')
        mult = math.factorial(a) * math.factorial(b_) * math.factorial(c) * math.factorial(d)
        shw = letter_shuffles(a, b_, c, d)
        const = {(0,) * e: Q2(math.factorial(e))}
        const = shuffle_sym(const, {(1,) * f: Q2(math.factorial(f))})
        for n, co in cu.items():
            coef = Q2(co.a, co.b) if hasattr(co, 'a') else Q2(Fr(str(co[0])), Fr(str(co[1])))
            acc = {}
            for w, m in shw.items():
                if not w:
                    sadd(acc, (), Q2(Fr(m * mult, n + 1)))
                    continue
                w1 = LETTERS[w[0]]
                fac = (ONE - w1 ** (n + 1)) / Q2(n + 1)
                if fac.iszero():
                    continue
                for word, cw in symG_words(w).items():
                    sadd(acc, word, cw * fac * Q2(m * mult))
            acc = shuffle_sym(acc, const)
            for word, cw in acc.items():
                sadd(blocks[cls], word, cw * coef * Q2(960))
    tot = {}
    for k, v in blocks.items():
        for word, c in v.items():
            sadd(tot, word, c)
        print('block %-9s pairing (order as built): %s   (reversed order): %s' % (k, pairing(v), pairing(v, True)))
    print('TOTAL W5/u pairing (order as built): %s   (reversed): %s' % (pairing(tot), pairing(tot, True)))
    print('paper: (-1485+1035 sqrt2, 900-540 sqrt2)')
    # depth-two correction C of the paper: Li_{3,2}(x,1) = G(0,0,1/x,0,1/x;1)
    def li32(x):
        ix = ONE / x
        key = ('0', '0', 'X', '0', 'X')
        seq = (ZERO,) + tuple(reversed((ZERO, ZERO, ix, ZERO, ix))) + (ONE,)
        return expand_words(symI(seq))
    bb = Q2(-1, 1)
    Cs = {}
    for x, cf in ((bb, Q2(-3600, 2160)), (Q2(0, Fr(1, 2)), Q2(11880, -8280))):
        for wd, c in li32(x).items():
            sadd(Cs, wd, c * cf)
    print('paper correction C = (-3600+2160 sqrt2) Li32(b,1) + (11880-8280 sqrt2) Li32(1/sqrt2,1): pairing %s (reversed %s)' % (pairing(Cs), pairing(Cs, True)))
