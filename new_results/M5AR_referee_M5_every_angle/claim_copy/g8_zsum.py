"""g8_zsum.py (copy of the notes' zsum.py, library only; see sec_g8_M2alpha) -- exact algebra of alternating nested sums (Z-sums), for symbolic-n identities.

Z_{a_1,...,a_k}(M) = sum_{M >= i_1 > i_2 > ... > i_k >= 1} prod_j sgn(a_j)^{i_j} / i_j^{|a_j|},   a_j in Z \\ {0},
Z_{()}(M) = 1.
An "expression" is a dict {word (tuple of nonzero ints): coefficient in Q(n)} meaning sum coef(n) Z_word(2n).
Coefficients live in the rational function field K = Q(n) (sympy sparse fields, exact).

Operations:
  stuffle(u, v)        : Z_u(M) Z_v(M) = sum Z_w(M)  (quasi-shuffle with a o b = sgn(a)sgn(b)(|a|+|b|), '+')
  shift(word, j)       : Z_word(2n + j) as an expression in Z_w(2n)   (j any integer; valid for 2n + j >= 0
                         and n large enough that no denominator 2n+t vanishes, t in the range used)
  harm_n(r) / harm_2n_minus_1(r) : H^(r)_{n-1}, H^(r)_{2n-1} as expressions
  hpl_coeff_terms(w)   : h_w(N) = sum coef eps^N / N^k Z_v(N-1)  (exact recursion from the HPL definition)
  hpl_coeff_at(w, delta): h_w(2n + delta) as an expression
"""
from fractions import Fraction as Fr
from functools import lru_cache
from sympy import QQ
from sympy.polys.fields import field

K, nn = field('n', QQ)
ONE = K(1)
ZERO = K(0)


def sgn(a):
    return 1 if a > 0 else -1


def circ(a, b):
    return sgn(a)*sgn(b)*(abs(a) + abs(b))


def add_into(E, word, cf):
    if cf == 0:
        return
    v = E.get(word)
    v = cf if v is None else v + cf
    if v == 0:
        E.pop(word, None)
    else:
        E[word] = v


def lin(*pairs):
    out = {}
    for cf, E in pairs:
        for w, v in E.items():
            add_into(out, w, cf*v)
    return out


@lru_cache(maxsize=None)
def stuffle(u, v):
    """dict {word: int multiplicity} for Z_u * Z_v (strict / Z-sum quasi-shuffle)."""
    if not u:
        return {v: 1}
    if not v:
        return {u: 1}
    out = {}
    a, up = u[0], u[1:]
    b, vp = v[0], v[1:]
    for w, m in stuffle(up, v).items():
        out[(a,) + w] = out.get((a,) + w, 0) + m
    for w, m in stuffle(u, vp).items():
        out[(b,) + w] = out.get((b,) + w, 0) + m
    for w, m in stuffle(up, vp).items():
        out[(circ(a, b),) + w] = out.get((circ(a, b),) + w, 0) + m
    return {w: m for w, m in out.items() if m != 0}


def mul(E1, E2):
    out = {}
    for w1, c1 in E1.items():
        for w2, c2 in E2.items():
            cc = c1*c2
            for w, m in stuffle(w1, w2).items():
                add_into(out, w, cc*m)
    return out


@lru_cache(maxsize=None)
def _shift(word, j):
    """Z_word(2n + j) as tuple of (word, coef) pairs."""
    if not word or j == 0:
        return ((word, ONE),)
    a, rest = word[0], word[1:]
    out = {word: ONE}
    sa = sgn(a); ea = abs(a)
    if j > 0:
        # Z_{a,rest}(2n+j) = Z_{a,rest}(2n) + sum_{t=1}^{j} sa^t/(2n+t)^ea Z_rest(2n+t-1)
        for t in range(1, j + 1):
            cf = K(sa**t)/(2*nn + t)**ea
            for w, c in _shift(rest, t - 1):
                add_into(out, w, cf*c)
    else:
        # Z_{a,rest}(2n+j) = Z_{a,rest}(2n) - sum_{t=j+1}^{0} sa^t/(2n+t)^ea Z_rest(2n+t-1)
        for t in range(j + 1, 1):
            cf = -K(sa**(t % 2))/(2*nn + t)**ea
            for w, c in _shift(rest, t - 1):
                add_into(out, w, cf*c)
    return tuple(out.items())


def shift(word, j):
    return dict(_shift(tuple(word), j))


def harm_2n(r):
    """H^(r)_{2n} = Z_r(2n)."""
    return {(r,): ONE}


def harm_2n_minus_1(r):
    """H^(r)_{2n-1} = Z_r(2n) - 1/(2n)^r."""
    return {(r,): ONE, (): -ONE/(2*nn)**r}


def harm_n_minus_1(r):
    """H^(r)_{n-1} = 2^{r-1}(Z_r(2n) + Z_{-r}(2n)) - 1/n^r."""
    return {(r,): K(2**(r - 1)), (-r,): K(2**(r - 1)), (): -ONE/nn**r}


@lru_cache(maxsize=None)
def hpl_coeff_terms(w):
    """h_w(N) = sum coef * eps^N / N^k * Z_v(N-1); returns tuple of (coef(Fraction), eps, k, v)."""
    w = tuple(w)
    if len(w) == 1:
        a = w[0]
        assert a in (1, -1)
        # h_a(N) = a^{N-1}/N = a * a^N / N
        return ((Fr(a), a, 1, ()),)
    inner = hpl_coeff_terms(w[1:])
    a = w[0]
    out = []
    if a == 0:
        for cf, eps, k, v in inner:
            out.append((cf, eps, k + 1, v))
    else:
        # h_{a,w}(N) = (1/N) sum_{i<=N-1} a^{N-1-i} h_w(i) = a * a^N/N * Z_{(a eps) k, v}(N-1)
        for cf, eps, k, v in inner:
            out.append((cf*a, a, 1, (a*eps*k,) + v))
    return tuple(out)


def hpl_coeff_at(w, delta):
    """h_w(2n + delta) as an expression in Z_v(2n) with Q(n) coefficients."""
    out = {}
    for cf, eps, k, v in hpl_coeff_terms(tuple(w)):
        pref = K(cf.numerator)/K(cf.denominator)*K(eps**(delta % 2))/(2*nn + delta)**k
        for ww, c in _shift(v, delta - 1):
            add_into(out, ww, pref*c)
    return out


# ------------------------------------------------------------ numeric (exact) evaluation, for testing
@lru_cache(maxsize=None)
def Z_val(word, M):
    if not word:
        return Fr(1)
    if M <= 0:
        return Fr(0)
    a, rest = word[0], word[1:]
    # Z_{a,rest}(M) = Z_{a,rest}(M-1) + sgn(a)^M / M^|a| Z_rest(M-1)
    return Z_val(word, M - 1) + Fr(sgn(a)**M, M**abs(a))*Z_val(rest, M - 1)


