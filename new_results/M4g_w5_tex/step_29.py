chk = r'''#!/usr/bin/env python
"""check_supplement.py:
(a) parse every coefficient and argument BACK out of supplement_W5.tex and
    compare exactly with the JSON (every term must round-trip);
(b) evaluate W5 from the parsed-back data at 400+ digits (direct series for
    Li_k and Li_{3,2}; Li5(-1) = -15 zeta(5)/16; L(chi8,5) via Hurwitz zeta)
    and compare with W5_5400.json (expect < 1e-390).
"""
import json, re
from fractions import Fraction
import mpmath as mp

BASE = '../M4e_w5_lchi'
terms = json.load(open(BASE + '/planner/R5_nice_reduced.json'))['terms']
n5json = json.load(open(BASE + '/cobracket.json'))['correction']
W5_STR = json.load(open(BASE + '/W5_5400.json'))['value']

ARG_FORMS = {
    'b': (Fraction(-1), Fraction(1)),
    '-b': (Fraction(1), Fraction(-1)),
    '-b^{2}': (Fraction(-3), Fraction(2)),
    '\\tfrac{1}{\\sqrt2}': (Fraction(0), Fraction(1, 2)),
    '\\tfrac12': (Fraction(1, 2), Fraction(0)),
    '1-\\tfrac{1}{\\sqrt2}': (Fraction(1), Fraction(-1, 2)),
    '2-\\sqrt2': (Fraction(2), Fraction(-1)),
    '2b': (Fraction(-2), Fraction(2)),
    '\\tfrac{2-\\sqrt2}{4}': (Fraction(1, 2), Fraction(-1, 4)),
    '\\tfrac{2+\\sqrt2}{4}': (Fraction(1, 2), Fraction(1, 4)),
    '-\\tfrac{b}{2}': (Fraction(1, 2), Fraction(-1, 2)),
    '\\tfrac{4-3\\sqrt2}{8}': (Fraction(1, 2), Fraction(-3, 8)),
    '1': (Fraction(1), Fraction(0)),
    '-1': (Fraction(-1), Fraction(0)),
}

def parse_rat(t):
    m = re.fullmatch(r'(-?)\\tfrac\{(\d+)\}\{(\d+)\}', t)
    if m:
        v = Fraction(int(m.group(2)), int(m.group(3)))
        return -v if m.group(1) == '-' else v
    return Fraction(int(t))

def parse_coeff(body, sign):
    if body.startswith('\\bigl('):
        assert body.endswith('\\bigr)')
        inner = body[6:-6]
        m = re.fullmatch(r'(.+?)([+-])(.+?)\\sqrt2', inner)
        p, q = parse_rat(m.group(1)), parse_rat(m.group(3))
        if m.group(2) == '-':
            q = -q
    elif body.endswith('\\sqrt2'):
        p, q = Fraction(0), parse_rat(body[:-6])
    else:
        p, q = parse_rat(body), Fraction(0)
    if sign == '-':
        p, q = -p, -q
    return (p, q)

def canon_field(pq):
    p, q = pq
    if q == 0:
        return str(p)
    if p == 0:
        return '%s*sqrt2' % q
    return '(%s%+s*sqrt2)' % (p, q)

def parse_atom(s):
    if s == 'L(\\chi_{8},5)':
        return 'Lchi8_5'
    m = re.fullmatch(r'\\operatorname\{Li\}_\{3,2\}\\bigl\((.+?),(.+?)\\bigr\)', s)
    if m:
        return ('li32', ARG_FORMS[m.group(1)], ARG_FORMS[m.group(2)])
    keys = []
    for part in s.split('\\,'):
        m = re.fullmatch(r'\\lambda(?:\^\{(\d+)\})?', part)
        if m:
            keys += ['lam'] * int(m.group(1) or 1); continue
        m = re.fullmatch(r'L(?:\^\{(\d+)\})?', part)
        if m:
            keys += ['L'] * int(m.group(1) or 1); continue
        if part == '\\pi^{2}':
            keys.append('pi2'); continue
        if part == '\\zeta(3)':
            keys.append('z3'); continue
        if part == '\\zeta(5)':
            keys.append('z5'); continue
        m = re.fullmatch(r'\\operatorname\{Li\}_\{(\d)\}\\bigl\((.+)\\bigr\)(?:\^\{(\d+)\})?', part)
        if m:
            k, a, n = m.group(1), m.group(2), int(m.group(3) or 1)
            name = 'Li%s(b)' % k if a == 'b' else 'Li%s(%s)' % (k, canon_field(ARG_FORMS[a]))
            keys += [name] * n; continue
        raise ValueError('unparsed atom part: %r' % part)
    return '*'.join(keys)

# ---------------- (a) round-trip -------------------------------------------
rt_terms = {}
rt_n5 = []
for line in open('supplement_W5.tex'):
    line = line.rstrip('\n')
    m = re.match(r'^(?:N_5=\{\})?&([+-])(.*?)\\\\$', line)
    if not m:
        continue
    sign, body = m.group(1), m.group(2)
    coeff_str, atom_str = body.split('\\,', 1)
    pq = parse_coeff(coeff_str, sign)
    atom = parse_atom(atom_str)
    if isinstance(atom, tuple):
        rt_n5.append({'x': [str(atom[1][0]), str(atom[1][1])],
                      'y': [str(atom[2][0]), str(atom[2][1])],
                      'coefficient': [str(pq[0]), str(pq[1])]})
    else:
        assert atom not in rt_terms, atom
        rt_terms[atom] = pq

ref_terms = {k: (Fraction(v[0]), Fraction(v[1])) for k, v in terms.items()}
assert rt_terms == ref_terms, set(rt_terms.items()) ^ set(ref_terms.items())
print('round-trip (a): all %d classical terms match the JSON exactly' % len(rt_terms))
ref_n5 = [{'x': [str(Fraction(e['x'][0])), str(Fraction(e['x'][1]))],
           'y': [str(Fraction(e['y'][0])), str(Fraction(e['y'][1]))],
           'coefficient': [str(Fraction(e['coefficient'][0])), str(Fraction(e['coefficient'][1]))]}
          for e in n5json]
assert rt_n5 == ref_n5, (rt_n5, ref_n5)
print('round-trip (a): all %d N5 terms match cobracket.json exactly' % len(rt_n5))

# ---------------- (b) evaluation at 400 digits ------------------------------
mp.dps = 430
sq2 = mp.sqrt(2)

def mpq(fr):
    return mp.mpf(fr.numerator) / fr.denominator

def fval(pq):
    return mpq(pq[0]) + mpq(pq[1]) * sq2

def li_series(k, x):
    N = int(450 / float(-mp.log10(abs(x)))) + 60
    s, xp = mp.mpf(0), x
    for j in range(1, N + 1):
        s += xp / mp.mpf(j) ** k
        xp *= x
    return s

def li32(x, y):
    N = int(450 / float(-mp.log10(abs(x)))) + 60
    s, xn, yj, h = mp.mpf(0), x, mp.mpf(1), mp.mpf(0)
    for n in range(2, N + 1):
        yj *= y
        h += yj / mp.mpf(n - 1) ** 2   # h = H_{n-1}(y)
        xn *= x
        s += xn / mp.mpf(n) ** 3 * h
    return s

LCHI = (mp.zeta(5, mp.mpf(1) / 8) + mp.zeta(5, mp.mpf(7) / 8)
        - mp.zeta(5, mp.mpf(3) / 8) - mp.zeta(5, mp.mpf(5) / 8)) / 8 ** 5

cache = {}
def atom_value(key):
    if key in cache:
        return cache[key]
    if key == 'Lchi8_5':
        v = LCHI
    else:
        v = None
        for f in split_key(key):
            m = re.fullmatch(r'Li(\d)\((.*)\)', f)
            base = {'lam': mp.log(2), 'L': mp.log(1 + sq2), 'pi2': mp.pi ** 2,
                    'z3': mp.zeta(3), 'z5': mp.zeta(5)}
            if f in base:
                fv = base[f]
            elif m:
                k = int(m.group(1))
                a = m.group(2)
                x = sq2 - 1 if a == 'b' else fval(parse_field_str(a))
                if x == -1:
                    fv = -mp.mpf(15) / 16 * mp.zeta(5)
                else:
                    fv = li_series(k, x)
            else:
                raise ValueError(f)
            v = fv if v is None else v * fv
    cache[key] = v
    return v

def split_key(key):
    parts, depth, cur = [], 0, ''
    for ch in key:
        if ch == '*' and depth == 0:
            parts.append(cur); cur = ''; continue
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        cur += ch
    parts.append(cur)
    return parts

def parse_field_str(s):
    s = s.replace('(', '').replace(')', '').strip()
    p, q = Fraction(0), Fraction(0)
    for t in re.findall(r'[+-]?[^+-]+', s):
        t = t.strip()
        if t.endswith('*sqrt2'):
            q += Fraction(t[:-6])
        else:
            p += Fraction(t)
    return (p, q)

W5est = mp.mpf(0)
for key, pq in rt_terms.items():
    W5est += (mpq(pq[0]) + mpq(pq[1]) * sq2) * atom_value(key)
print('classical part evaluated:', mp.nstr(W5est, 8))
for e in rt_n5:
    c = mpq(Fraction(e['coefficient'][0])) + mpq(Fraction(e['coefficient'][1])) * sq2
    x = mpq(Fraction(e['x'][0])) + mpq(Fraction(e['x'][1])) * sq2
    y = mpq(Fraction(e['y'][0])) + mpq(Fraction(e['y'][1])) * sq2
    W5est += c * li32(x, y)
target = mp.mpf(W5_STR)
diff = W5est - target
print('W5 estimate :', mp.nstr(W5est, 20))
print('W5 (5400bit):', mp.nstr(target, 20))
print('difference  :', mp.nstr(diff, 5))
assert abs(diff) < mp.mpf('1e-390'), abs(diff)
print('CHECK (b): |W5est - W5| < 1e-390  PASS (got ~%s)' % mp.nstr(abs(diff), 3))
'''
open('check_supplement.py','w').write(chk)
import subprocess
r=subprocess.run(['nice','-n','10','<project>/.venv/bin/python','check_supplement.py'],capture_output=True,text=True,timeout=1500)
print(r.stdout)
print(r.stderr[-1500:])
