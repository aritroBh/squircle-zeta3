
reflib_algebra = '''"""reflib.py -- HEXC2R referee's INDEPENDENT library. No imports from HEXC2/HEXC1 code.

Level-6 CMZV word algebra, written from the definitions:
  letters: -1 = dt/t ; j in 0..5 -> omega_j = zeta^j dt/(1-zeta^j t), zeta = e^{i pi/3}.
  word w = (a_1,...,a_n), FIRST LETTER OUTERMOST:
  I(w; 0->x) = int_{x>t_1>...>t_n>0} omega_{a_1}(t_1)...omega_{a_n}(t_n).
  admissible (convergent as written): first letter != 0, last letter != -1.
"""
from fractions import Fraction as Fr
from functools import cache
from itertools import product
from math import gcd
from functools import reduce

N = 6
LETTERS = tuple(range(-1, N))

def vadd(*terms):
    """terms are (dict, coef); exact rational vector add."""
    out = {}
    for d, c in terms:
        for w, a in d.items():
            nw = out.get(w, Fr(0)) + c * a
            if nw: out[w] = nw
            elif w in out: del out[w]
    return out

@cache
def sh(u, v):
    """shuffle product of two words -> dict."""
    if not u: return {v: Fr(1)}
    if not v: return {u: Fr(1)}
    out = {}
    for w, a in sh(u[1:], v).items():
        w2 = (u[0],) + w
        out[w2] = out.get(w2, Fr(0)) + a
    for w, a in sh(u, v[1:]).items():
        w2 = (v[0],) + w
        out[w2] = out.get(w2, Fr(0)) + a
    return out

def mul(U, V):
    """shuffle product of two vectors."""
    out = {}
    for u, a in U.items():
        for v, b in V.items():
            for w, c in sh(u, v).items():
                nw = out.get(w, Fr(0)) + a * b * c
                if nw: out[w] = nw
                elif w in out: del out[w]
    return out

@cache
def reg(w):
    """shuffle-regularized finite part of the word value, T=0 convention:
    strip leading 0's and trailing -1's via the shuffle recursion."""
    w = tuple(w)
    if not w: return {(): Fr(1)}
    if w[0] != 0 and w[-1] != -1: return {w: Fr(1)}
    if w[0] == 0:
        single, tail = (0,), w[1:]
    else:
        single, tail = (-1,), w[:-1]
    terms = sh(single, tail)
    c = terms.pop(w)
    return vadd(*[(reg(t), Fr(-a, c)) for t, a in terms.items()])

def regularize(D):
    return vadd(*[(reg(w), a) for w, a in D.items()])

@cache
def words(n, admissible=True):
    if n == 0: return ((),)
    return tuple(w for w in product(LETTERS, repeat=n)
                 if w[-1] != -1 and (not admissible or w[0] != 0))

# ---- series (stuffle) side ----
@cache
def to_series(w):
    """word -> ((s_1,c_1),...,(s_r,c_r)): dt/t increments s; letter j closes a block
    (s, (j - prev_j) mod 6) with prev_0 = 0. Fails if the word ends with -1."""
    out = []; s = 1; last = 0
    for j in w:
        if j == -1: s += 1
        else:
            out.append((s, (j - last) % N)); last = j; s = 1
    if s != 1: raise ValueError("word ends with -1")
    return tuple(out)

@cache
def from_series(srs):
    out = []; last = 0
    for s, c in srs:
        last = (last + c) % N
        out.extend([-1] * (s - 1) + [last])
    return tuple(out)

@cache
def stuffle(u, v):
    """stuffle product of series-index tuples."""
    if not u: return {v: Fr(1)}
    if not v: return {u: Fr(1)}
    out = {}
    for w, a in stuffle(u[1:], v).items():
        w2 = (u[0],) + w
        out[w2] = out.get(w2, Fr(0)) + a
    for w, a in stuffle(u, v[1:]).items():
        w2 = (v[0],) + w
        out[w2] = out.get(w2, Fr(0)) + a
    merged = (u[0][0] + v[0][0], (u[0][1] + v[0][1]) % N)
    for w, a in stuffle(u[1:], v[1:]).items():
        w2 = (merged,) + w
        out[w2] = out.get(w2, Fr(0)) + a
    return out

def ds(u, v):
    """regularized double-shuffle relation row: regularize(stuffle(u,v) - shuffle(u,v)).
    u may start with 0 (the series side still converges); v admissible."""
    a = {}
    for w, c in stuffle(to_series(u), to_series(v)).items():
        a[from_series(w)] = c
    return regularize(vadd((a, 1), (sh(u, v), -1)))

def distribution(n, d):
    """Zhao's regularized distribution relation of order d | 6, for words over
    letters {-1} u {0, d, 2d, ...} (residues multiples of d)."""
    rows = []
    # powers of the letter-sum C = sum_{j in (N/d)Z/NZ...} used for the leading-0 correction
    C = {(j,): Fr(1) for j in range(N // d, N, N // d)}
    cp = [{(): Fr(1)}]
    for k in range(1, n + 1):
        cp.append(vadd((mul(cp[-1], C), Fr(1, k))))
    for w in product((-1,) + tuple(range(0, N, d)), repeat=n):
        w = tuple(w)
        if w[-1] == -1: continue
        # lift: each letter j (multiple of d) -> one of {j/d + k*N/d : k=0..d-1}; -1 -> -1
        lift_choices = []
        for j in w:
            if j == -1:
                lift_choices.append(((-1,),))
            else:
                lift_choices.append(tuple(((j // d + k * (N // d)) % N,) for k in range(d)))
        left = {}
        for tup in product(*lift_choices):
            w2 = tuple(x[0] for x in tup)
            left[w2] = left.get(w2, Fr(0)) + d ** w.count(-1)
        right = reg(w)
        for k in range(1, n + 1):
            if w[k - 1] != 0: break
            right = vadd((right, 1), (mul(cp[k], reg(w[k:])), 1))
        rows.append(regularize(vadd((left, 1), (right, -1))))
    return rows

def weight1_rows():
    """branch relations among Li_1(zeta^j) = -Log(1-zeta^j) (principal branch):
    (1)+(5) = 0  [ (1-z)(1-z^5) = 1, args cancel ]
    (1)-(5)-2(2)+2(4) = 0  [ arg combination ]."""
    r1 = vadd(({(1,): Fr(1)}, 1), ({(5,): Fr(1)}, 1))
    r2 = vadd(({(1,): Fr(1)}, 1), ({(5,): Fr(1)}, -1), ({(2,): Fr(1)}, -2), ({(4,): Fr(1)}, 2))
    return [r1, r2]

def primitive(row):
    row = {w: Fr(c) for w, c in row.items() if c}
    if not row: return ()
    den = reduce(lambda a, b: a * b // gcd(a, b), (c.denominator for c in row.values()), 1)
    vals = sorted((w, c.numerator * (den // c.denominator)) for w, c in row.items())
    g = reduce(gcd, (c for _, c in vals))
    if vals[0][1] < 0: g = -g
    return tuple((w, c // g) for w, c in vals)

def make_relations(n, lower):
    """all proved relation rows of weight n (as primitive integer tuples), with family counts.
    lower: dict k -> rows of weight k (for lifts)."""
    rows = []; seen = set(); counts = {}
    def emit(row, fam):
        r = primitive(row)
        if not r: return
        counts[fam + '_nonzero'] = counts.get(fam + '_nonzero', 0) + 1
        if r not in seen:
            seen.add(r); rows.append(r)
            counts[fam + '_new'] = counts.get(fam + '_new', 0) + 1
    if n == 1:
        for r in weight1_rows():
            rr = primitive(r)
            if rr and rr not in seen:
                seen.add(rr); rows.append(rr)
                counts['weight1_branch_new'] = counts.get('weight1_branch_new', 0) + 1
    else:
        for a in range(1, n):
            b = n - a
            for u in words(a, False):
                u_adm = u[0] != 0
                if u_adm and a > b: continue
                for v in words(b):
                    if u_adm and a == b and u > v: continue
                    emit(ds(u, v), 'double_shuffle')
    for d in (2, 3, 6):
        for row in distribution(n, d): emit(row, 'distribution')
    for k, basis in lower.items():
        for row in basis:
            for w in words(n - k):
                emit(mul(dict(row), {w: Fr(1)}), 'lift')
    return rows, counts
'''
with open("reflib.py", "w") as f:
    f.write(reflib_algebra)
print("written", len(reflib_algebra), "chars")
