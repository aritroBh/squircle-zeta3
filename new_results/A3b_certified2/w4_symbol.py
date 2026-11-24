"""Exact, constant-retaining symbol check for main.tex, Conjecture conj:W4.

Uses the supplied QS2 arithmetic; does not call its floating-point S-unit finder.
G_c(w;e) means integration from c to e, with the first letter outermost.
Coincident differences use unit tangents (log 0 has regularized value zero).
The check is at symbol level only: pi^2 and zeta(3) have zero full symbol.
"""
import sys
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import product
from pathlib import Path

PAPER = Path('<project>/paper/squircle_note')
sys.dont_write_bytecode = True
sys.path.insert(0, str(PAPER / 'scripts'))
from qsqrt2 import QS2 as Q, SQ2 as R, A, B

ZERO, ONE = Q(0), Q(1)


@cache
def logvec(z):
    """log|z| = j/2 log2 + m log a, found and verified without floats."""
    z = abs(Q.c(z))
    if z == ZERO:
        return (F(0), F(0))  # unit-tangent constant term, not an ordinary log
    norm = abs(z.norm())
    p, q = norm.numerator, norm.denominator
    if p & (p - 1) or q & (q - 1):
        raise ValueError(f'not an S-unit: {z}, norm={norm}')
    j = p.bit_length() - q.bit_length()
    u = z / R**j
    m = 0
    while u > ONE:
        u /= A
        m += 1
    while u < ONE:
        u *= A
        m -= 1
        if u > ONE:
            raise ValueError(f'not an S-unit: {z}')
    assert z == R**j * A**m
    return (F(j, 2), F(m))


def clean(d):
    return {k: v for k, v in d.items() if v != ZERO}


def add(*ts):
    out = defaultdict(lambda: ZERO)
    for t in ts:
        for k, v in t.items():
            out[k] += v
    return clean(out)


def scale(t, c):
    return clean({k: v * c for k, v in t.items()})


@cache
def shuffles(u, v):
    if not u:
        return (v,)
    if not v:
        return (u,)
    return tuple((u[0],) + w for w in shuffles(u[1:], v)) + tuple(
        (v[0],) + w for w in shuffles(u, v[1:]))


def shuffle(s, t):
    out = defaultdict(lambda: ZERO)
    for u, c in s.items():
        for v, d in t.items():
            for w in shuffles(u, v):
                out[w] += c * d
    return clean(out)


def log_symbol(z):
    return {(i,): Q(v) for i, v in enumerate(logvec(z)) if v}


def concat(s, t):
    out = defaultdict(lambda: ZERO)
    for u, c in s.items():
        for v, d in t.items():
            out[u + v] += c * d
    return clean(out)


@cache
def g_symbol(c, w, e):
    if not w:
        return {(): ONE}
    if c == e:
        return {}
    out = {}
    points = (e,) + w + (c,)
    for i, x in enumerate(w, 1):
        ell = add(log_symbol(x - points[i - 1]),
                  scale(log_symbol(x - points[i + 1]), -1))
        out = add(out, concat(g_symbol(c, w[:i - 1] + w[i:], e), ell))
    return out


def li_symbol(n, z):
    out = scale(log_symbol(ONE - z), -1)
    for _ in range(n - 1):
        out = concat(out, log_symbol(z))
    return out


# l(v) = G_1(b;v) + G_1(-a;v) - G_1(0;v) - log2;
# m(v) = G_1(a;v) + G_1(-b;v) - G_1(0;v) - log2.
# Polynomial keys (word, k) denote G_1(word;v) * log(2)^k.
LT = {((B,), 0): ONE, ((-A,), 0): ONE, ((ZERO,), 0): -ONE,
      ((), 1): -ONE}
MT = {((A,), 0): ONE, ((-B,), 0): ONE, ((ZERO,), 0): -ONE,
      ((), 1): -ONE}


def log_product(*factors):
    out = {((), 0): ONE}
    for f in factors:
        nxt = defaultdict(lambda: ZERO)
        for (u, k), c in out.items():
            for (v, j), d in f.items():
                for w in shuffles(u, v):
                    nxt[(w, k + j)] += c * d
        out = clean(nxt)
    return out


def w4_g_expression():
    """Exact finite GPL expression, endpoints e=a minus e=b; no quadrature."""
    out = defaultdict(lambda: ZERO)
    for z, factors, coefficient in [
            (ZERO, (LT, LT, LT), (8 - 7 * R) / 48),
            (-A, (LT, LT, LT), Q(F(-1, 3))),
            (ZERO, (LT, LT, MT), R / 16)]:
        for (w, k), d in log_product(*factors).items():
            out[((z,) + w, k)] += coefficient * d
    return clean(out)


def w4_symbol():
    out = {}
    for (w, k), c in w4_g_expression().items():
        s = add(g_symbol(ONE, w, A), scale(g_symbol(ONE, w, B), -1))
        for _ in range(k):
            s = shuffle(s, log_symbol(Q(2)))
        out = add(out, scale(s, c))
    return out


def q(p=0, r=0):
    return Q(F(p), F(r))


def conjecture_terms():
    """Literal transcription of conj:W4. (coefficient, tuple of atoms)."""
    l, h, p, z = ('log', A), ('log', Q(2)), ('pi2',), ('z3',)
    d = ('li', 2, B)
    t, u, v = ('li', 3, B), ('li', 3, -B), ('li', 3, R / 2)
    terms = []
    def term(c, *atoms):
        terms.append((c, atoms))
    for arg, c in [(B, q(1, -573)), (-B, q(-1, -573)), (-B*B, q(0, 56)),
                   (R/2, q(0, 17)), (q('1/2'), q('-37/4', '-17/16')),
                   (ONE-R/2, q(10, 86)), (2-R, q(-10, 86))]:
        term(c, ('li', 4, arg))
    for c, a in [(q(40, 7), t), (q(-40, 7), u), (q(-72), v),
                 (q('63/8', '-3/16'), z)]:
        term(c, l, a)
    for c, a in [(q(6, 13), t), (q(6, -13), u), (q(0, '-47/2'), v),
                 (q('-3/8', '329/128'), z)]:
        term(c, h, a)
    term(q(30), d, d)
    term(q(30, '7/2'), d, l, l)
    term(q(-24), d, h, l)
    term(q(0, '-19/8'), d, h, h)
    term(q('-19/4', '1/24'), d, p)
    term(q('103/384', '-4963/5760'), p, p)
    term(q('67/24', '185/96'), p, l, l)
    term(q('-35/12', '41/12'), p, h, l)
    term(q('5/12', '-1105/384'), p, h, h)
    term(q('13/24', '-25/24'), l, l, l, l)
    term(q('-73/6', '-29/3'), l, l, l, h)
    term(q('63/8', '67/16'), l, l, h, h)
    term(q('41/24', '17/6'), l, h, h, h)
    term(q('-61/192', '137/192'), h, h, h, h)
    return terms


def atom_symbol(a):
    if a[0] == 'li':
        return li_symbol(a[1], a[2])
    if a[0] == 'log':
        return log_symbol(a[1])
    assert a in [('pi2',), ('z3',)]
    return {}


def rhs_symbol():
    out = {}
    for c, atoms in conjecture_terms():
        s = {(): ONE}
        for atom in atoms:
            s = shuffle(s, atom_symbol(atom))
        out = add(out, scale(s, c))
    return out


def check_tex_transcription():
    # Reuse only the supplied parser functions; importing its script runs numerics.
    import ast
    import re
    import sympy as sp
    from sympy.parsing.sympy_parser import (parse_expr,standard_transformations,
        implicit_multiplication_application,convert_xor)
    tree = ast.parse((PAPER/'scripts'/'check_tex_formulas.py').read_text())
    tree.body = [node for node in tree.body if isinstance(node,ast.FunctionDef)
                 and node.name in {'group','token','tex2py'}]
    scope = {'re':re}
    exec(compile(tree,'supplied_tex_parser','exec'),scope)
    tex = (PAPER/'main.tex').read_text().split(r'\label{conj:W4}',1)[1]
    tex = tex.split(r'\begin{multline*}',1)[1].split(r'\end{multline*}',1)[0]
    Ls,Hs,Zs,bs = sp.symbols('L lam z3 b')
    li = {i:sp.Function(f'Li{i}') for i in (2,3,4)}
    parsed = parse_expr(scope['tex2py'](tex.replace('W_4=','')),
        local_dict={'L':Ls,'lam':Hs,'z3':Zs,'b':bs,'sqrt':sp.sqrt,'pi':sp.pi,
                    **{f'Li{i}':li[i] for i in li}},
        transformations=standard_transformations+(implicit_multiplication_application,convert_xor))
    def atom(a):
        if a[0] == 'li':
            return li[a[1]](a[2].sympy())
        if a[0] == 'log':
            x,y = logvec(a[1])
            return sp.Rational(x.numerator,x.denominator)*Hs+sp.Rational(y.numerator,y.denominator)*Ls
        return sp.pi**2 if a == ('pi2',) else Zs
    target = sum(c.sympy()*sp.prod(atom(a) for a in atoms) for c,atoms in conjecture_terms())
    assert sp.expand(parsed.subs(bs,B.sympy())-target) == 0


def ban_float():
    def forbidden(self):
        raise AssertionError('a proof checker tried converting QS2 to float')
    Q.__float__ = forbidden


def self_check():
    # Each elementary GPL test uses a different construction on the RHS.
    for n in range(1, 5):
        for z in [B, -B, -B*B, R/2, Q(F(1, 2)), ONE-R/2, 2-R]:
            got = scale(g_symbol(ZERO, (ZERO,) * (n-1) + (ONE/z,), ONE), -1)
            assert got == li_symbol(n, z), (n, z, got)
    # Repeated letters, path reversal, endpoint regularization, and shuffle.
    for c, e in [(ZERO, ONE), (ONE, A), (ONE, B)]:
        for x in [ZERO, A, -A, B, -B]:
            s = g_symbol(c, (x,), e)
            for n in range(2, 5):
                s = shuffle(s, g_symbol(c, (x,), e))
                import math
                assert s == scale(g_symbol(c, (x,) * n, e), math.factorial(n))
        for x, y in product([ZERO, A, -A, B, -B], repeat=2):
            assert shuffle(g_symbol(c, (x,), e), g_symbol(c, (y,), e)) == add(
                g_symbol(c, (x, y), e), g_symbol(c, (y, x), e))


def main():
    ban_float()
    self_check()
    check_tex_transcription()
    lhs, rhs = w4_symbol(), rhs_symbol()
    diff = add(lhs, scale(rhs, -1))
    print('Exact GPL/Li_n/shuffle/unit-tangent checks: PASS')
    print('Exact comparison to the formula parsed from main.tex: PASS')
    print('GPL terms (each evaluated at a minus b):', len(w4_g_expression()))
    print('Symbol alphabet: log2, log(1+sqrt2)')
    print('word, W4 coefficient, conjecture coefficient')
    for w in product(range(2), repeat=4):
        print(''.join('2' if x == 0 else 'a' for x in w), lhs.get(w, ZERO), rhs.get(w, ZERO))
    print('Difference:', diff)
    assert not diff
    print('PROVED-MODULO-FULL-SYMBOL-KERNEL: S(W4) = S(conjecture RHS)')
    if '--export-gpl' in sys.argv:
        lines = ['# Sum c*log(2)^k*(G_1(word;a)-G_1(word;b)); first letter outermost.',
                 'coefficient_p\tcoefficient_q\tk\tword']
        for (w,k),c in w4_g_expression().items():
            lines.append(f'{c.p}\t{c.q}\t{k}\t{w}')
        Path('gpl_expression.tsv').write_text('\n'.join(lines)+'\n')


if __name__ == '__main__':
    main()
