import shutil, subprocess, os
shutil.copy('supplement_W5.tex','supplement_W5_old.tex')
gen = r'''#!/usr/bin/env python
"""gen_supplement.py [TERMS_JSON [OUT_TEX]]: build the W5 supplement table from a terms JSON.

argv[1]: terms JSON with key 'terms' mapping atom name -> [p, q]
         (default ../M4e_w5_lchi/planner/R5_nice_reduced.json)
argv[2]: output .tex path (default supplement_W5.tex)
N5 is read from ../M4e_w5_lchi/cobracket.json; the W5/M4 digit strings from
W5_5400.json and M4_800.txt. No coefficient is hand-typed.
"""
import json, re, math, sys
from fractions import Fraction

BASE = '../M4e_w5_lchi'
TERMS_PATH = sys.argv[1] if len(sys.argv) > 1 else BASE + '/planner/R5_nice_reduced.json'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'supplement_W5.tex'
terms = json.load(open(TERMS_PATH))['terms']
n5 = json.load(open(BASE + '/cobracket.json'))['correction']
W5 = json.load(open(BASE + '/W5_5400.json'))['value']
M4 = open('../M4a_decomposition/M4_800.txt').read().strip()

def parse_field(s):
    """e.g. '(-1+1*sqrt2)' -> (Fraction(-1), Fraction(1)); tolerant of stray parens."""
    s = s.replace('(', '').replace(')', '').strip()
    p, q = Fraction(0), Fraction(0)
    for t in re.findall(r'[+-]?[^+-]+', s):
        t = t.strip()
        if t.endswith('*sqrt2'):
            q += Fraction(t[:-6])
        else:
            p += Fraction(t)
    return (p, q)

# Closed forms of every argument used, as exact Q(sqrt2) pairs.
ARG_FORMS = [   # (latex, (p, q)) with arg = p + q*sqrt(2)
    ('b',                     (Fraction(-1), Fraction(1))),
    ('-b',                    (Fraction(1), Fraction(-1))),
    ('-b^{2}',                (Fraction(-3), Fraction(2))),
    ('b^{2}',                 (Fraction(3), Fraction(-2))),
    ('\\tfrac{1}{\\sqrt2}',   (Fraction(0), Fraction(1, 2))),
    ('-\\tfrac{1}{\\sqrt2}',  (Fraction(0), Fraction(-1, 2))),
    ('\\tfrac12',             (Fraction(1, 2), Fraction(0))),
    ('1-\\tfrac{1}{\\sqrt2}', (Fraction(1), Fraction(-1, 2))),
    ('2-\\sqrt2',             (Fraction(2), Fraction(-1))),
    ('2b',                    (Fraction(-2), Fraction(2))),
    ('\\tfrac{2-\\sqrt2}{4}', (Fraction(1, 2), Fraction(-1, 4))),
    ('\\tfrac{2+\\sqrt2}{4}', (Fraction(1, 2), Fraction(1, 4))),
    ('-\\tfrac{b}{2}',        (Fraction(1, 2), Fraction(-1, 2))),
    ('\\tfrac{4-3\\sqrt2}{8}',(Fraction(1, 2), Fraction(-3, 8))),
    ('b^{4}',                 (Fraction(17), Fraction(-12))),
    ('1-b^{4}',               (Fraction(-16), Fraction(12))),
    ('1',                     (Fraction(1), Fraction(0))),
    ('-1',                    (Fraction(-1), Fraction(0))),
]
ARG_LATEX = {pq: lx for lx, pq in ARG_FORMS}
SQ2 = math.sqrt(2.0)

def arg_latex(pq):
    assert pq in ARG_LATEX, 'no closed form registered for %r' % (pq,)
    return ARG_LATEX[pq]

def arg_factor(a):
    if a == 'b':
        return 'b'
    return arg_latex(parse_field(a))

# numeric cross-check: every JSON field element against its registered closed form
for key in terms:
    for m in re.finditer(r'Li(\d)\(((?:[^()]|\([^()]*\))*)\)', key):
        if m.group(2) != 'b':
            pq = parse_field(m.group(2))
            print('arg check  %-22s -> %-20s float %.17g'
                  % (m.group(2), arg_latex(pq), float(pq[0]) + float(pq[1]) * SQ2))
for e in n5:
    for xy in ('x', 'y'):
        pq = (Fraction(e[xy][0]), Fraction(e[xy][1]))
        print('N5 %s check [%s,%s] -> %-8s float %.17g'
              % (xy, e[xy][0], e[xy][1], arg_latex(pq), float(pq[0]) + float(pq[1]) * SQ2))

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

def factor_latex(f, n):
    base = {'lam': '\\lambda', 'L': 'L', 'pi2': '\\pi^{2}',
            'z3': '\\zeta(3)', 'z5': '\\zeta(5)'}
    if f in base:
        b = base[f]
        if f in ('lam', 'L') and n > 1:
            b += '^{%d}' % n
        if f == 'pi2' and n > 1:
            b = '\\pi^{%d}' % (2 * n)
        return b
    m = re.fullmatch(r'Li(\d)\((.*)\)', f)
    assert m, f
    out = '\\operatorname{Li}_{%s}\\bigl(%s\\bigr)' % (m.group(1), arg_factor(m.group(2)))
    if n > 1:
        out += '^{%d}' % n
    return out

def split_factors(key):
    parts, depth, cur = [], 0, ''
    for ch in key:
        if ch == '*':
            if depth == 0:
                parts.append(cur); cur = ''; continue
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        cur += ch
    parts.append(cur)
    return parts

def atom_latex(key):
    if key == 'Lchi8_5':
        return 'L(\\chi_{8},5)'
    parts, out, i = split_factors(key), [], 0
    while i < len(parts):
        n = 1
        while i + n < len(parts) and parts[i + n] == parts[i]:
            n += 1
        out.append(factor_latex(parts[i], n))
        i += n
    return '\\,'.join(out)

# ---------- group orders, built from the keys present -----------------------
W4_ARGS = [(-1, 1), (1, -1), (-3, 2), (0, Fraction(1, 2)), (Fraction(1, 2), 0),
           (1, Fraction(-1, 2)), (2, -1)]
LI5_EXTRA_ARGS = [(3, -2), (0, Fraction(-1, 2)), (-2, 2),
                  (Fraction(1, 2), Fraction(-1, 4)), (Fraction(1, 2), Fraction(1, 4)),
                  (Fraction(1, 2), Fraction(-1, 2)), (Fraction(1, 2), Fraction(-3, 8)),
                  (17, -12), (-16, 12)]
LI5_ARG_ORDER = [tuple(Fraction(x) for x in a) for a in W4_ARGS + LI5_EXTRA_ARGS]

def key_arg(key):
    m = re.fullmatch(r'Li(\d)\(((?:[^()]|\([^()]*\))*)\)', key)
    return parse_field(m.group(2))

li5_keys = [k for k in terms if k.startswith('Li5(')]
pos = {}
for k in li5_keys:
    pq = key_arg(k)
    assert pq in LI5_ARG_ORDER, 'Li5 argument not in the ordering list: %r' % k
    pos[k] = LI5_ARG_ORDER.index(pq)
LI5_ORDER = sorted(li5_keys, key=lambda k: pos[k])
LI5_ORDER += [k for k in ('z5', 'Lchi8_5') if k in terms]

LI4_ALL = [w + '*Li4(' + a + ')' for w in ('lam', 'L') for a in
           ('(-1+1*sqrt2)', '(1-1*sqrt2)', '(-3+2*sqrt2)', '1/2*sqrt2',
            '1/2', '(1-1/2*sqrt2)', '(2-1*sqrt2)')]
LI4_ORDER = [k for k in LI4_ALL if k in terms]
W3W2_ALL = [w + '*z3' for w in ('lam*lam', 'lam*L', 'L*L', 'pi2', 'Li2(b)')] + \
    [w + '*Li3(' + a + ')' for w in ('lam*lam', 'lam*L', 'L*L', 'pi2', 'Li2(b)')
     for a in ('(-1+1*sqrt2)', '(1-1*sqrt2)', '1/2*sqrt2')] + \
    [w + '*Li2(b)' for w in ('lam*lam*lam', 'lam*lam*L', 'lam*L*L', 'L*L*L',
                             'lam*pi2', 'L*pi2')] + \
    ['lam*Li2(b)*Li2(b)', 'L*Li2(b)*Li2(b)']
W3W2_ORDER = [k for k in W3W2_ALL if k in terms]
PURE_ALL = ['lam*lam*lam*lam*lam', 'lam*lam*lam*lam*L', 'lam*lam*lam*L*L',
            'lam*lam*L*L*L', 'lam*L*L*L*L', 'L*L*L*L*L',
            'lam*lam*lam*pi2', 'lam*lam*L*pi2', 'lam*L*L*pi2', 'L*L*L*pi2',
            'lam*pi2*pi2', 'L*pi2*pi2']
PURE_ORDER = [k for k in PURE_ALL if k in terms]
placed = LI5_ORDER + LI4_ORDER + W3W2_ORDER + PURE_ORDER
assert len(placed) == len(set(placed)) == len(terms) and set(placed) == set(terms), \
    set(terms) ^ set(placed)
print('groups: Li5 %d, Li4 %d, w3w2 %d, pure %d (total %d)'
      % (len(LI5_ORDER), len(LI4_ORDER), len(W3W2_ORDER), len(PURE_ORDER), len(placed)))

def term_line(key):
    p, q = Fraction(terms[key][0]), Fraction(terms[key][1])
    s, body = coeff_latex(p, q)
    return '&%s%s\\,%s' % (s, body, atom_latex(key))

def n5_line(e):
    p, q = Fraction(e['coefficient'][0]), Fraction(e['coefficient'][1])
    s, body = coeff_latex(p, q)
    xl = arg_latex((Fraction(e['x'][0]), Fraction(e['x'][1])))
    yl = arg_latex((Fraction(e['y'][0]), Fraction(e['y'][1])))
    return '&%s%s\\,\\operatorname{Li}_{3,2}\\bigl(%s,%s\\bigr)' % (s, body, xl, yl)

def fifty(v):
    ip, fp = v.strip().split('.')
    return ip + '.' + fp[:50]

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
def emit_align(lines):
    lines[0] = lines[0].replace('&+', '&', 1)   # no leading + on the first line
    A(r'\begin{align*}')
    A('\\\\\n'.join(lines))
    A(r'\end{align*}')
emit_align(['N_5={}' + n5_line(n5[0])] + [n5_line(e) for e in n5[1:]])
A(r'\section*{The classical part}')
A(r'The difference $W_5-N_5$ is the following $\mathbb{Q}(\sqrt2)$-combination. '
  r'First the $\operatorname{Li}_5$ terms, $\zeta(5)$ and $L(\chi_8,5)$:')
emit_align(['W_5-N_5={}' + term_line(LI5_ORDER[0])] + [term_line(k) for k in LI5_ORDER[1:]])
A(r'Then the products of a tetralogarithm with $\lambda$ or $L$:')
emit_align([term_line(k) for k in LI4_ORDER])
A(r'Then the products of weight three and weight two:')
emit_align([term_line(k) for k in W3W2_ORDER])
A(r'Finally the pure products in $\lambda$, $L$ and $\pi^2$:')
emit_align([term_line(k) for k in PURE_ORDER])
A(r'\section*{Status}')
A(r'Numerically')
A(r'\begin{gather*}')
A('W_5=%s\\ldots\\\\' % fifty(W5))
A('M_4=%s\\ldots' % fifty(M4))
A(r'\end{gather*}')
A(r'(50 digits each.) The identity displayed above is a numerical identity: the right-hand side agrees '
  r'with $W_5$ to 1596 digits, and its full point symbol and $(3,2)$ cobracket agree with those of $W_5$ '
  r'exactly. It is not a proof. The coefficient of $L(\chi_8,5)$ depends on the representative, because '
  r'$L(\chi_8,5)$ itself lies in the $\operatorname{Li}_5$ span.')
A(r'\end{document}')
open(OUT, 'w').write('\n'.join(L) + '\n')
print('wrote %s: %d lines' % (OUT, len(L)))
'''
open('gen_supplement.py','w').write(gen)
# check_supplement.py: argv + new ARG_FORMS + optional sign
c=open('check_supplement.py').read()
c=c.replace("import json, re\nfrom fractions import Fraction\nimport mpmath as mp",
            "import json, re, sys\nfrom fractions import Fraction\nimport mpmath as mp")
c=c.replace("terms = json.load(open(BASE + '/planner/R5_nice_reduced.json'))['terms']",
            "TERMS_PATH = sys.argv[1] if len(sys.argv) > 1 else BASE + '/planner/R5_nice_reduced.json'\n"
            "TEX_PATH = sys.argv[2] if len(sys.argv) > 2 else 'supplement_W5.tex'\n"
            "terms = json.load(open(TERMS_PATH))['terms']")
c=c.replace("    '-b^{2}': (Fraction(-3), Fraction(2)),",
            "    '-b^{2}': (Fraction(-3), Fraction(2)),\n    'b^{2}': (Fraction(3), Fraction(-2)),")
c=c.replace("    '\\\\tfrac{1}{\\\\sqrt2}': (Fraction(0), Fraction(1, 2)),",
            "    '\\\\tfrac{1}{\\\\sqrt2}': (Fraction(0), Fraction(1, 2)),\n    '-\\\\tfrac{1}{\\\\sqrt2}': (Fraction(0), Fraction(-1, 2)),")
c=c.replace("    '\\\\tfrac{4-3\\\\sqrt2}{8}': (Fraction(1, 2), Fraction(-3, 8)),",
            "    '\\\\tfrac{4-3\\\\sqrt2}{8}': (Fraction(1, 2), Fraction(-3, 8)),\n    'b^{4}': (Fraction(17), Fraction(-12)),\n    '1-b^{4}': (Fraction(-16), Fraction(12)),")
c=c.replace("m = re.match(r'^(?:N_5=\\{\\}|W_5-N_5=\\{\\})?&([+-])(.*?)(?:\\\\\\\\)?$', line)",
            "m = re.match(r'^(?:N_5=\\{\\}|W_5-N_5=\\{\\})?&([+-]?)(.*?)(?:\\\\\\\\)?$', line)\n    if m and m.group(2) == '':\n        m = (m.group(1), '+', m.group(3))")
# adjust tuple unpacking after the regex
c=c.replace("    sign, body = m.group(1), m.group(2)",
            "    sign, body = (m[0], m[1]) if isinstance(m, tuple) else (m.group(1), m.group(2))")
c=c.replace("for line in open('supplement_W5.tex'):", "for line in open(TEX_PATH):")
open('check_supplement.py','w').write(c)
print('files written')
