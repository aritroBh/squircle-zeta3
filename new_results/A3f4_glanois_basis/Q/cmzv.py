"""Exact shuffle algebra and Cayley pullbacks; H letters -1=dt/t, j=ζ^j dt/(1-ζ^j t)."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from functools import cache
from itertools import product, permutations
import json, csv, time
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / 'paper/squircle_note/scripts'))
sys.path.insert(0, str(HERE.parent.parent / 'A3_W4_proof'))
import w4_symbol as source
from qsqrt2 import QS2 as Q, A, B, SQ2


def log(s):
    s = time.strftime('%Y-%m-%d %H:%M:%S') + ' ' + s
    print(s, flush=True)
    with (HERE / 'progress.log').open('a') as f:
        f.write(s + '\n')


def add(*args):
    out = Counter()
    for d, c in args:
        for w, a in d.items(): out[w] += c*a
    return {w:a for w,a in out.items() if a}


@cache
def sh(u, v):
    if not u: return {v:1}
    if not v: return {u:1}
    out = Counter()
    for w,a in sh(u[1:],v).items(): out[(u[0],)+w] += a
    for w,a in sh(u,v[1:]).items(): out[(v[0],)+w] += a
    return dict(out)


def mul(u, v):
    out = Counter()
    for w,a in u.items():
        for z,b in v.items():
            for t,c in sh(w,z).items(): out[t] += a*b*c
    return {w:a for w,a in out.items() if a}


@cache
def reg(w):
    """Shuffle finite part H(0;1)=H(-1;1)=0, including trailing dt/t."""
    if not w: return {(): F(1)}
    if w[0] != 0 and w[-1] != -1: return {w:F(1)}
    if w[0] == 0: single, tail = (0,), w[1:]
    else: single, tail = (-1,), w[:-1]
    terms = dict(sh(single,tail))
    c = terms.pop(w)
    return add(*[(reg(t),F(-a,c)) for t,a in terms.items()])


def regularize(d):
    return add(*[(reg(w),a) for w,a in d.items()])


# Cayley exponents of the configuration, including infinity (None).
POINTS = {None:0, Q(0):4, Q(1):6, Q(-1):2, A:7, B:5, -A:1, -B:3}


@cache
def path(w,c,e):
    """G_c(w;e), root exponents, via 0; scale each radial path to endpoint 1."""
    if c == e: return {} if w else {():F(1)}
    out = {}
    for k in range(len(w)+1):
        left = tuple((e-x)%8 for x in w[:k])
        right = tuple((c-x)%8 for x in reversed(w[k:]))
        out = add((out,1),(sh(left,right),(-1)**k))
    return regularize(out)


def forms_path(forms,c,e):
    out = {}
    for choices in product(*forms):
        w = tuple(x for x,a in choices)
        coefficient = 1
        for x,a in choices: coefficient *= a
        out = add((out,1),(path(w,c,e),coefficient))
    return out


def dlog(p, q=None):
    # Pullback dlog((x-p)/(x-q)); infinity maps to root exponent 0.
    return ((POINTS[p],1),(POINTS[q],-1))


def det(x,y):
    if x is None: return Q(0 if y is None else 1)
    if y is None: return Q(-1)
    return x-y


def cross(x,p0,p1,pinf):
    denominator = det(x,pinf)*det(p1,p0)
    return None if denominator == Q(0) else det(x,p0)*det(p1,pinf)/denominator


@cache
def find_cross(z):
    for p0,p1,pinf,e in permutations(POINTS,4):
        if p0 is None or e is None: continue
        if pinf is not None and min(p0,e) < pinf < max(p0,e): continue
        if cross(e,p0,p1,pinf) == z:
            assert cross(p0,p0,p1,pinf) == Q(0)
            return p0,p1,pinf,e
    raise ValueError(('no cross ratio',z))


@cache
def atom(a):
    if a[0] == 'log':
        if a[1] == Q(2): return {(4,):F(-1)}
        assert a[1] == A
        return forms_path([dlog(Q(0))],POINTS[Q(1)],POINTS[A])
    if a[0] == 'pi2':
        ipi = {(2,):F(2),(6,):F(-2)}
        return add((mul(ipi,ipi),-1))
    if a[0] == 'z3': return {(-1,-1,0):F(1)}
    _,n,z = a
    p0,p1,pinf,e = find_cross(z)
    return add((forms_path([dlog(p0,pinf)]*(n-1)+[dlog(p1,pinf)],
                          POINTS[p0],POINTS[e]),-1))


def vectors():
    original = source.w4_g_expression()
    with (HERE.parent.parent/'A3_W4_proof/gpl_expression.tsv').open() as f:
        rows = list(csv.DictReader((line for line in f if not line.startswith('#')),delimiter='\t'))
    exported = {(repr(w),k):(c.p,c.q) for (w,k),c in original.items()}
    assert len(rows) == len(original)
    for row in rows:
        assert exported[row['word'],int(row['k'])] == (F(row['coefficient_p']),F(row['coefficient_q']))
    log(f'GPL input rows verified against source: {len(rows)}')
    lhs = [{},{}]
    for (w,k),coefficient in original.items():
        forms = [dlog(p) for p in w]
        term = add((forms_path(forms,POINTS[Q(1)],POINTS[A]),1),
                   (forms_path(forms,POINTS[Q(1)],POINTS[B]),-1))
        for _ in range(k): term = mul(term,atom(('log',Q(2))))
        for j,c in enumerate((coefficient.p,coefficient.q)):
            lhs[j] = add((lhs[j],1),(term,c))
    rhs = [{},{}]
    for coefficient,atoms in source.conjecture_terms():
        term = {():F(1)}
        for a in atoms: term = mul(term,atom(a))
        for j,c in enumerate((coefficient.p,coefficient.q)):
            rhs[j] = add((rhs[j],1),(term,c))
    difference = [add((lhs[j],1),(rhs[j],-1)) for j in range(2)]
    for name,vec in [('W',lhs),('R',rhs),('D',difference)]:
        log(f'{name} supports (rational, sqrt2): {tuple(map(len,vec))}')
    data = {name:[[[list(w),str(c)] for w,c in sorted(v.items())] for v in vec]
            for name,vec in [('W',lhs),('R',rhs),('D',difference)]}
    (HERE/'vectors.json').write_text(json.dumps(data))
    return lhs,rhs,difference


if __name__ == '__main__':
    log('Exact Cayley conversion started; single process; no source writes.')
    for a in dict.fromkeys(a for c,ats in source.conjecture_terms() for a in ats):
        if a[0] == 'li': log(f'Li cross ratio (p0,p1,pinf,e), n={a[1]}, z={a[2]}: {find_cross(a[2])}')
    vectors()
