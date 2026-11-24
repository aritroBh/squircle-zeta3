gen = r'''#!/usr/bin/env python
"""gen_supplement.py: build supplement_W5.tex ONLY from the JSON inputs.

Inputs (read-only):
  ../M4e_w5_lchi/planner/R5_nice_reduced.json  (key 'terms')
  ../M4e_w5_lchi/cobracket.json                (key 'correction', the N5 terms)
  ../M4e_w5_lchi/W5_5400.json                  (key 'value', W5 to ~1600 digits)
  ../M4a_decomposition/M4_800.txt              (M4 to 800 digits)
No coefficient is hand-typed; every displayed number comes from the JSON.
"""
import json, re, sys
from fractions import Fraction

BASE = '../M4e_w5_lchi'
terms = json.load(open(BASE + '/planner/R5_nice_reduced.json'))['terms']
n5 = json.load(open(BASE + '/cobracket.json'))['correction']
W5 = json.load(open(BASE + '/W5_5400.json'))['value']
M4 = open('../M4a_decomposition/M4_800.txt').read().strip()

# ---------- field elements p + q*sqrt2 -------------------------------------
def parse_field(s):
    """'(-1+1*sqrt2)' -> (Fraction(-1), Fraction(1)); also '1/2*sqrt2', '1/2'."""
    s = s.strip()
    if s.startswith('(') and s.endswith(')'):
        s = s[1:-1]
    p, q = Fraction(0), Fraction(0)
    for t in re.findall(r'[+-]?[^+-]+', s):
        t = t.strip()
        if t.endswith('*sqrt2'):
            q += Fraction(t[:-6])
        else:
            p += Fraction(t)
    return (p, q)

# The 12 Li5 / 7 Li4 / 3 Li3 arguments, as closed forms (exact Q(sqrt2) pairs).
ARG_FORMS = [   # (latex, (p, q)) with arg = p + q*sqrt(2)
    ('b',                    (Fraction(-1), Fraction(1))),
    ('-b',                   (Fraction(1), Fraction(-1))),
    ('-b^{2}',               (Fraction(-3), Fraction(2))),
    ('\\tfrac{1}{\\sqrt2}',  (Fraction(0), Fraction(1, 2))),
    ('\\tfrac12',            (Fraction(1, 2), Fraction(0))),
    ('1-\\tfrac{1}{\\sqrt2}',(Fraction(1), Fraction(-1, 2))),
    ('2-\\sqrt2',            (Fraction(2), Fraction(-1))),
    ('2b',                   (Fraction(-2), Fraction(2))),
    ('\\tfrac{2-\\sqrt2}{4}',(Fraction(1, 2), Fraction(-1, 4))),
    ('\\tfrac{2+\\sqrt2}{4}',(Fraction(1, 2), Fraction(1, 4))),
    ('-\\tfrac{b}{2}',       (Fraction(1, 2), Fraction(-1, 2))),
    ('\\tfrac{4-3\\sqrt2}{8}',(Fraction(1, 2), Fraction(-3, 8))),
    ('1',                    (Fraction(1), Fraction(0))),
    ('-1',                   (Fraction(-1), Fraction(0))),
]
ARG_LATEX = {pq: lx for lx, pq in ARG_FORMS}

def arg_latex(pq):
    assert pq in ARG_LATEX, 'no closed form registered for %r' % (pq,)
    return ARG_LATEX[pq]

# numeric cross-check of each registered form against the JSON field elements
import math
SQ2 = math.sqrt(2.0)
def field_float(pq):
    return float(pq[0]) + float(pq[1]) * SQ2
for key in list(terms):
    for m in re.finditer(r'Li\d\(([^)]*)\)', key):
        a = m.group(1)
        if a == 'b':
            continue
        pq = parse_field(a)
        lx = arg_latex(pq)
        print('arg check %-24s -> %-22s float %.17g' % (a, lx, field_float(pq)))
for e in n5:
    for xy in ('x', 'y'):
        pq = parse_field('(%s%+s*sqrt2)' % (xy == 'x' and e['x'][0] or e['y'][0],
                                            0)) if False else (Fraction(e[xy][0]), Fraction(e[xy][1]))
        print('N5 %s check [%s,%s] -> %s float %.17g' % (xy, e[xy][0], e[xy][1], arg_latex(pq), field_float(pq)))

# ---------- coefficients ----------------------------------------------------
def fmt_rat(r):
    assert r > 0
    if r.denominator == 1:
        return str(r.numerator)
    return '\\tfrac{%d}{%d}' % (r.numerator, r.denominator)

def coeff_latex(p, q):
    """(p,q) -> (sign, body); overall sign factored out, paper style."""
    s = '+'
    if p < 0 or (p == 0 and q < 0):
        s, p, q = '-', -p, -q
    if q == 0:
        return s, fmt_rat(p)
    if p == 0:
        return s, fmt_rat(abs(q)) + '\\sqrt2'
    inner = fmt_rat(p) + ('+' if q > 0 else '-') + fmt_rat(abs(q)) + '\\sqrt2'
    return s, '\\bigl(%s\\bigr)' % inner

# ---------- atoms -----------------------------------------------------------
def factor_latex(f, n):
    base = {'lam': '\\lambda', 'L': 'L', 'pi2': '\\pi^{2}',
            'z3': '\\zeta(3)', 'z5': '\\zeta(5)'}
    if f in base:
        b = base[f]
        if f in ('lam', 'L') and n > 1:
            b += '^{%d}' % n
        return b
    m = re.fullmatch(r'Li(\d)\((.*)\)', f)
    assert m, f
    a = m.group(2)
    al = 'b' if a == 'b' else arg_latex(parse_field(a))
    return '\\operatorname{Li}_{%s}\\bigl(%s\\bigr)' % (m.group(1), al)

def atom_latex(key):
    if key == 'Lchi8_5':
        return 'L(\\chi_{8},5)'
    parts, out, i = key.split('*'), [], 0
    while i < len(parts):
        n = 1
        while i + n < len(parts) and parts[i + n] == parts[i]:
            n += 1
        out.append(factor_latex(parts[i], n))
        i += n
    return '\\,'.join(out)

# ---------- groups and their order ------------------------------------------
LI5_ORDER = ['Li5((-1+1*sqrt2))', 'Li5((1-1*sqrt2))', 'Li5((-3+2*sqrt2))',
             'Li5(1/2*sqrt2)', 'Li5(1/2)', 'Li5((1-1/2*sqrt2))',
             'Li5((2-1*sqrt2))', 'Li5((-2+2*sqrt2))', 'Li5((1/2-1/4*sqrt2))',
             'Li5((1/2+1/4*sqrt2))', 'Li5((1/2-1/2*sqrt2))',
             'Li5((1/2-3/8*sqrt2))', 'z5', 'Lchi8_5']
LI4_ORDER = [w + '*Li4(' + a + ')' for w in ('lam', 'L') for a in
             ('(-1+1*sqrt2)', '(1-1*sqrt2)', '(-3+2*sqrt2)', '1/2*sqrt2',
              '1/2', '(1-1/2*sqrt2)', '(2-1*sqrt2)')]
W3W2_ORDER = [w + '*z3' for w in ('lam*lam', 'lam*L', 'L*L', 'pi2', 'Li2(b)')] + \
    [w + '*Li3(' + a + ')' for w in ('lam*lam', 'lam*L', 'L*L', 'pi2', 'Li2(b)')
     for a in ('(-1+1*sqrt2)', '(1-1*sqrt2)', '1/2*sqrt2')] + \
    [w + '*Li2(b)' for w in ('lam*lam*lam', 'lam*lam*L', 'lam*L*L', 'L*L*L',
                             'lam*pi2', 'L*pi2')] + \
    ['lam*Li2(b)*Li2(b)', 'L*Li2(b)*Li2(b)']
PURE_ORDER = ['lam*lam*lam*lam*lam', 'lam*lam*lam*lam*L', 'lam*lam*lam*L*L',
              'lam*lam*L*L*L', 'lam*L*L*L*L', 'L*L*L*L*L',
              'lam*lam*lam*pi2', 'lam*lam*L*pi2', 'lam*L*L*pi2', 'L*L*L*pi2',
              'lam*pi2*pi2', 'L*pi2*pi2']
assert set(LI5_ORDER + LI4_ORDER + W3W2_ORDER + PURE_ORDER) == set(terms), \
    set(terms) ^ set(LI5_ORDER + LI4_ORDER + W3W2_ORDER + PURE_ORDER)
assert len(terms) == 68 and len(LI5_ORDER) == 14 and len(LI4_ORDER) == 14 \
    and len(W3W2_ORDER) == 28 and len(PURE_ORDER) == 12

def term_line(key):
    p, q = Fraction(terms[key][0]), Fraction(terms[key][1])
    s, body = coeff_latex(p, q)
    return '&%s%s\\,%s\\\\' % (s, body, atom_latex(key))

def n5_line(e):
    p, q = Fraction(e['coefficient'][0]), Fraction(e['coefficient'][1])
    s, body = coeff_latex(p, q)
    xl = arg_latex((Fraction(e['x'][0]), Fraction(e['x'][1])))
    yl = arg_latex((Fraction(e['y'][0]), Fraction(e['y'][1])))
    return '&%s%s\\,\\operatorname{Li}_{3,2}\\bigl(%s,%s\\bigr)\\\\' % (s, body, xl, yl)

def fifty(v):
    v = v.strip()
    ip, fp = v.split('.')
    return ip + '.' + fp[:50]

# ---------- the document ----------------------------------------------------
L = []
A = L.append
A(r'\documentclass[11pt]{article}')
A(r'\usepackage{amsmath,amssymb}')
A(r'\usepackage[margin=1in]{geometry}')
A(r'\allowdisplaybreaks')
A(r'\title{Supplement: the weight-five part $W_5$ of the squircle coefficient $M_4$}')
A(r'\author{}')
A(r'\date{}')
A(r'\begin{document}')
A(r'\maketitle')
A(r'\section*{Notation}')
A(r'Throughout, $b=\sqrt2-1$, $\lambda=\log2$ and $L=\log(1+\sqrt2)$. The double polylogarithm is '
  r'$\operatorname{Li}_{s,t}(x,y)=\sum_{n>m>0}x^ny^m/(n^sm^t)$, and $\operatorname{Li}_n(z)=\sum_{k\ge1}z^k/k^n$. '
  r'$\chi_8$ is the real character modulo $8$, with $\chi_8(1)=\chi_8(7)=1$, $\chi_8(3)=\chi_8(5)=-1$ and '
  r'$\chi_8(n)=0$ for $n$ even, and $L(\chi_8,5)=\sum_{n\ge1}\chi_8(n)/n^5$. '
  r'Every coefficient below is an element $p+q\sqrt2$ of $\mathbb{Q}(\sqrt2)$ with $p,q$ rational.')
A(r'\section*{The depth-two part $N_5$}')
A(r'The $(3,2)$ cobracket of $W_5$ fixes four double-polylogarithm terms exactly:')
A(r'\begin{align*}')
A('N_5={}' + n5_line(n5[0]))
for e in n5[1:]:
    A(n5_line(e))
A(r'\end{align*}')
A(r'\section*{The classical part}')
A(r'The difference $W_5-N_5$ is the following $\mathbb{Q}(\sqrt2)$-combination. '
  r'First the $\operatorname{Li}_5$ terms, $\zeta(5)$ and $L(\chi_8,5)$:')
A(r'\begin{align*}')
for k in LI5_ORDER:
    A(term_line(k))
A(r'\end{align*}')
A(r'Then the products of a tetralogarithm with $\lambda$ or $L$:')
A(r'\begin{align*}')
for k in LI4_ORDER:
    A(term_line(k))
A(r'\end{align*}')
A(r'Then the products of weight three and weight two:')
A(r'\begin{align*}')
for k in W3W2_ORDER:
    A(term_line(k))
A(r'\end{align*}')
A(r'Finally the pure products in $\lambda$, $L$ and $\pi^2$:')
A(r'\begin{align*}')
for k in PURE_ORDER:
    A(term_line(k))
A(r'\end{align*}')
A(r'\section*{Status}')
A(r'Numerically $W_5=%s\ldots$ and $M_4=%s\ldots$ (50 digits each).' % (fifty(W5), fifty(M4)))
A(r'The identity displayed above is a numerical identity: the right-hand side agrees with $W_5$ to 1596 digits, '
  r'and its full point symbol and $(3,2)$ cobracket agree with those of $W_5$ exactly. It is not a proof. '
  r'The coefficient of $L(\chi_8,5)$ depends on the representative, because $L(\chi_8,5)$ itself lies in the '
  r'$\operatorname{Li}_5$ span.')
A(r'\end{document}')
open('supplement_W5.tex', 'w').write('\n'.join(L) + '\n')
print('wrote supplement_W5.tex: %d lines' % len(L))
'''
open('gen_supplement.py','w').write(gen)
print('gen_supplement.py written', len(gen))
