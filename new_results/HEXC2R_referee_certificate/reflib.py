"""reflib.py -- HEXC2R referee's INDEPENDENT library. No imports from HEXC2/HEXC1 code.

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
    terms = dict(sh(single, tail))  # copy: sh is cached, pop would mutate the cache
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


# ================= independent numeric evaluator =================
# Values from the DEFINITION as iterated sums. For a word with letters = forms
# (dt/t or a*dt/(1-a*t)), scanned into blocks ((s_1,a_1),...,(s_r,a_r)):
#     I(w; 0->x) = sum_{n_1>...>n_r>=1}  x^{n_1} prod_i a_i^{n_i-n_{i+1}} / prod_i n_i^{s_i}
# (n_{r+1}=0). With |a_i| <= 1 and |x| = 1/2 every summand is O(2^{-n_1} poly), so the
# sum is evaluated by truncation. Equivalently with z_1 = a_1 x, z_i = a_i/a_{i-1}:
#     value = sum_{n_2>=1} z_2^{n_2}/n_2^{s_2} * L1(n_2) * Q(n_2),
#     L1(n) = sum_{n_1>n} z_1^{n_1}/n_1^{s_1},  Q = inner multiple sum (blocks 3..r),
# all accumulated incrementally in O(K) operations, K ~ (dps+20)/log10(2).
import mpmath as mp

_STATE2 = {}

def _blocks(letters):
    """letters: tuple of ('T',) or ('A', alpha). Returns block tuple ((s,alpha),...)."""
    out = []; s = 1
    for L in letters:
        if L[0] == 'T':
            s += 1
        else:
            out.append((s, L[1])); s = 1
    if s != 1:
        raise ValueError("divergent word (trailing dt/t)")
    return tuple(out)

def _zpow_list(z, K):
    out = [mp.mpc(1)] * (K + 1)
    for i in range(1, K + 1): out[i] = out[i - 1] * z
    return out

def _inv_tables(K):
    """inv[n] = (n^{-1}, n^{-2}, n^{-3}, n^{-4}) for n = 1..K, cached per precision."""
    key = ('inv', K, mp.mp.dps)
    if key in _STATE2: return _STATE2[key]
    tab = [None] * (K + 1)
    for n in range(1, K + 1):
        m = mp.mpf(n)
        tab[n] = (1/m, 1/m**2, 1/m**3, 1/m**4)
    _STATE2[key] = tab
    return tab

def msum(blocks, x):
    """the truncated multiple sum, |a_i|<=1, |x|<=1/2. Geometric tail 2^{-K}."""
    r = len(blocks)
    if r == 0: return mp.mpc(1)
    dps = mp.mp.dps
    K = int((dps + 20) / 0.30103) + 30
    inv = _inv_tables(K)
    s = [b[0] for b in blocks]
    a = [mp.mpc(b[1]) for b in blocks]
    x = mp.mpc(x)
    z1 = a[0] * x
    z1p = _zpow_list(z1, K)
    i1 = [mp.mpf(1)] + [inv[n][s[0] - 1] for n in range(1, K + 1)]
    L = [mp.mpc(0)] * (K + 2)
    for n1 in range(K, 0, -1):
        L[n1 - 1] = L[n1] + z1p[n1] * i1[n1]
    if r == 1:
        return L[0]
    z2 = a[1] / a[0]
    z2p = _zpow_list(z2, K)
    i2 = [mp.mpf(1)] + [inv[n][s[1] - 1] for n in range(1, K + 1)]
    total = mp.mpc(0)
    if r == 2:
        for n2 in range(1, K):
            total += z2p[n2] * i2[n2] * L[n2]
        return total
    if r == 3:
        z3 = a[2] / a[1]
        z3p = _zpow_list(z3, K)
        i3 = [mp.mpf(1)] + [inv[n][s[2] - 1] for n in range(1, K + 1)]
        Q = mp.mpc(0)
        for n2 in range(1, K):
            total += z2p[n2] * i2[n2] * L[n2] * Q
            Q += z3p[n2] * i3[n2]
        return total
    if r == 4:
        z3 = a[2] / a[1]; z4 = a[3] / a[2]
        z3p = _zpow_list(z3, K); z4p = _zpow_list(z4, K)
        i3 = [mp.mpf(1)] + [inv[n][s[2] - 1] for n in range(1, K + 1)]
        i4 = [mp.mpf(1)] + [inv[n][s[3] - 1] for n in range(1, K + 1)]
        P = mp.mpc(0); Q = mp.mpc(0)
        for n2 in range(1, K):
            total += z2p[n2] * i2[n2] * L[n2] * Q
            Q += z3p[n2] * i3[n2] * P
            P += z4p[n2] * i4[n2]
        return total
    raise ValueError("depth > 4 not needed here")

# ---- CMZV words (letters -1, 0..5) ----
def _zeta():
    return mp.e ** (mp.pi * mp.j / 3)

def _aforms():
    Z = _zeta()
    return [None] + [ -Z**j / (1 - Z**j) for j in range(1, 6) ]

def word_letters_A(w):
    """CMZV word -> form letters for the 0 -> x side."""
    Z = _zeta()
    return tuple(('T',) if j == -1 else ('A', Z**j) for j in w)

def word_letters_B(u):
    """CMZV word u (for the x -> 1 side) -> (sign, form letters) after s = 1-t reflection.
    Derivation: I(u; 1/2->1) = I(ref(u); 0->1/2) with ref = REVERSED word, forms g(s)=f(1-s):
    letter -1: g = ds/(1-s) = (A,1) ; letter 0: g = ds/s = ('T') ;
    letter j: g = -a_j ds/(1-a_j s), a_j = -zeta^j/(1-zeta^j)  (per-letter sign -1).
    Overall sign = (-1)^{#{letters 1..5}}.  (No (-1)^|u| prefactor: the dt -> -ds signs
    are absorbed by the flip of each iterated limit.)"""
    AF = _aforms()
    sign = 1
    out = []
    for j in reversed(u):
        if j == -1:
            out.append(('A', mp.mpc(1)))
        elif j == 0:
            out.append(('T',))
        else:
            out.append(('A', AF[j])); sign = -sign
    return sign, tuple(out)

_CV_CACHE = {}
def cmzv_value(w):
    """regularized value of an ADMISSIBLE level-6 CMZV word (first letter outermost)."""
    w = tuple(w)
    if not w: return mp.mpc(1)
    assert w[0] != 0 and w[-1] != -1, "word not admissible"
    key = (w, mp.mp.dps)
    if key in _CV_CACHE: return _CV_CACHE[key]
    tot = mp.mpc(0)
    for k in range(len(w) + 1):
        pre, suf = w[:k], w[k:]
        A = msum(_blocks(word_letters_A(suf)), mp.mpf(1) / 2)
        if pre:
            sgn, letters = word_letters_B(pre)
            B = sgn * msum(_blocks(letters), mp.mpf(1) / 2)
        else:
            B = mp.mpc(1)
        tot += B * A
    _CV_CACHE[key] = tot
    return tot

def cmzv_combo(D):
    """D: dict word -> Fraction coefficient."""
    tot = mp.mpc(0)
    for w, c in D.items():
        tot += mp.mpf(c.numerator) / c.denominator * cmzv_value(w)
    return tot

# ---- HPL words (letters 0, +/-1 = dt/t, dt/(1-t), dt/(1+t)), first letter outermost ----
def hpl_value(w, x, _memo={}):
    """Regularized HPL value H_w(x), letters 0,±1 (0 = dt/t, ±1 = dt/(1∓t)), first letter
    outermost. Trailing 0's (log x divergence) are removed by MY shuffle identity
    H_u * H_0 = sum_v mult_v H_v  (H_0(x) = log x, T = 0 convention); convergent words are
    evaluated as truncated multiple sums (geometric for |x| <= 1/2)."""
    w = tuple(w)
    key = (w, mp.mp.dps)
    if key in _memo: return _memo[key]
    if not w: return mp.mpc(1)
    L = mp.log(x)
    if w == (0,):
        return L
    if w[-1] != 0:
        # letter 1 -> +dt/(1-t) = (A,1); letter -1 -> +dt/(1+t) = -(A,-1): sign flip per -1
        letters = tuple(('T',) if j == 0 else ('A', mp.mpc(j)) for j in w)
        sgn = (-1) ** sum(1 for j in w if j == -1)
        val = sgn * msum(_blocks(letters), mp.mpf(x))
    else:
        u = w[:-1]
        shuf = dict(sh(u, (0,)))        # H_u ш H_(0) = sum_v mult_v H_v
        m = shuf.pop(w)                 # multiplicity of w = (trailing zero run)+1
        tot = L * hpl_value(u, x)
        for v, mult in shuf.items():
            tot -= mult * hpl_value(v, x)
        val = tot / m
    _memo[key] = val
    return val


# ================= independent conversion machinery (Task 2) =================
# My own derivation (see step_06 header). Notation: letter (j) = zeta^j du/(1-zeta^j u),
# dlog(u - zeta^k) = -letter((-k) mod 6). MAPS[mi] = (e0, eh, e1, em, einf) = exponents of
# f(0), f(1/2), f(1), f(-1), f(inf) for the 12 maps.
from itertools import product as _product

MAPS6 = [
    (0,1,2,5,4), (0,5,4,1,2), (1,2,3,0,5), (1,0,5,2,3),
    (2,1,0,3,4), (2,3,4,1,0), (3,2,1,4,5), (3,4,5,2,1),
    (4,5,0,3,2), (4,3,2,5,0), (5,0,1,4,3), (5,4,3,0,1),
]

def _dlogU(k):
    """dlog(u - zeta^k) as a letter vector: -letter((-k) mod 6)."""
    return {((-k) % 6,): Fr(-1)}

def hpl_form_pullback(mi, sig):
    """pullback of HPL form omega_sig (0: dlog x; 1: -dlog(x-1); -1: dlog(x+1)) under g=f^{-1},
    as a level-6 letter vector. g*dlog(x-q) = dlog(u-f(q)) - dlog(u-f(inf))."""
    e0, eh, e1, em, einf = MAPS6[mi]
    if sig == 0:
        return vadd((_dlogU(e0), 1), (_dlogU(einf), -1))
    if sig == 1:
        return vadd((_dlogU(e1), -1), (_dlogU(einf), 1))
    if sig == -1:
        return vadd((_dlogU(em), 1), (_dlogU(einf), -1))
    raise ValueError(sig)

def _shift6(w, s):
    return tuple((-1 if j == -1 else (j + s) % 6) for j in w)

def seg_value(w, c, e):
    """I(w; zeta^c -> zeta^e), first letter outermost, shuffle-regularized (T=0).
    My path composition through 0 (the outer letters w[:k] go with the e-side):
        I(w; a->b) = sum_k I(w[:k]; 0->b) * I(w[k:]; a->0),
        I(v; a->0) = (-1)^|v| I(reversed(v); 0->a),
    products of regularized values = shuffle products of the vectors."""
    w = tuple(w)
    if c == e:
        return {(): Fr(1)} if not w else {}
    out = {}
    n = len(w)
    for k in range(n + 1):
        left = regularize({_shift6(w[:k], e): Fr(1)})          # I(w[:k]; 0 -> zeta^e)
        right = regularize({_shift6(tuple(reversed(w[k:])), c): Fr(1)})  # I(rev(w[k:]); 0 -> zeta^c)
        out = vadd((out, 1), (mul(left, right), Fr((-1) ** (n - k))))
    return regularize(out)

def forms_value(formvecs, c, e):
    """iterated integral of a tuple of letter-vectors (form combinations) from zeta^c to zeta^e."""
    out = {}
    for choices in _product(*[tuple(fv.items()) for fv in formvecs]):
        w = tuple(x for cw, a in choices for x in cw)
        coef = Fr(1)
        for cw, a in choices: coef *= a
        out = vadd((out, 1), (seg_value(w, c, e), coef))
    return out

def my_hpl_to_cmzv(w, mi=0):
    """H_w(1/2) (letters 0,±1, first letter outermost) as a level-6 CMZV vector."""
    w = tuple(w)
    if not w: return {(): Fr(1)}
    formvecs = [hpl_form_pullback(mi, sig) for sig in w]
    e0, eh, e1, em, einf = MAPS6[mi]
    return forms_value(formvecs, e0, eh)

# ---- cross-ratio frames for Li atoms (my own enumeration) ----
def _crq(e, p0, p1, pinf):
    """rational cross-ratio (e,p0;p1,pinf); points are Fractions or None (=inf)."""
    def sub(a, b):
        if a is None: return None
        if b is None: return None
        return a - b
    def div(a, b):
        if a is None and b is None: return Fr(1)
        if a is None: return None
        if b is None: return Fr(0)
        return a / b
    f1 = div(sub(e, p0), sub(e, pinf)); f2 = div(sub(p1, pinf), sub(p1, p0))
    if f1 is None or f2 is None: return None
    return f1 * f2

def cross_reps(z):
    """all (e, p0, p1, pinf) with cr = z, from {0, 1/2, 1, -1, inf}."""
    from itertools import permutations as _perm
    PTS = [Fr(0), Fr(1, 2), Fr(1), Fr(-1), None]
    out = []
    for e, p0, p1, pinf in _perm(PTS, 4):
        if _crq(e, p0, p1, pinf) == Fr(z):
            out.append((e, p0, p1, pinf))
    return out

def my_li_atom(n, z, mi=0, rep=0, frame=None):
    """Li_n(z) as a level-6 CMZV vector under map mi, using the cross-ratio frame
    z = (e,p0;p1,pinf):  Li_n(z) = -G(0^{n-1},1; z) pulled back to u:
    outer forms dt/t -> dlog(u-f(p0)) - dlog(u-f(pinf));
    inner form dt/(t-1) -> dlog(u-f(p1)) - dlog(u-f(pinf)).  Path u: f(p0) -> f(e).
    frame: explicit (e,p0,p1,pinf) (with None = inf) to force a representation."""
    if frame is not None:
        e, p0, p1, pinf = frame
        assert _crq(e, p0, p1, pinf) == Fr(z), (frame, z)
    else:
        e, p0, p1, pinf = cross_reps(z)[rep]
    e0, eh, e1, em, einf = MAPS6[mi]
    EXP = {Fr(0): e0, Fr(1, 2): eh, Fr(1): e1, Fr(-1): em, None: einf}
    def dlx(q):
        # dlog(x - q) pulls back to dlog(u - f(q)) - dlog(u - f(inf)); q = inf gives 0-form
        if q is None:
            return {}
        return vadd((_dlogU(EXP[q]), 1), (_dlogU(EXP[None]), -1))
    f0 = vadd((dlx(p0), 1), (dlx(pinf), -1))
    # G-letter 1 = form dt/(t-1): pulls back to dlog(u-f(p1)) - dlog(u-f(pinf))
    f1 = vadd((dlx(p1), 1), (dlx(pinf), -1))
    formvecs = [f0] * (n - 1) + [f1]
    cexp = EXP[p0] if p0 is not None else einf
    eexp = EXP[e] if e is not None else einf
    return {w: -c for w, c in forms_value(formvecs, cexp, eexp).items()}

# ---- elementary atoms (my own derivations) ----
def my_atom_log2():  return {(3,): Fr(-1)}              # value (3,) = -log(1-zeta^3) = -log 2
def my_atom_log3():  return {(2,): Fr(-1), (4,): Fr(-1)}  # (2,)+(4,) = -log((1-z^2)(1-z^4)) = -log 3
def my_atom_log32(): return vadd((my_atom_log3(), 1), (my_atom_log2(), -1))
def my_atom_ipi():   return {(1,): Fr(3)}               # value (1,) = -Log(1-zeta) = i pi/3
def my_atom_zeta2(): return {(-1, 0): Fr(1)}
def my_atom_zeta3(): return {(-1, -1, 0): Fr(1)}
def my_atom_zeta4(): return {(-1, -1, -1, 0): Fr(1)}
def my_atom_pi2():   return mul(my_atom_ipi(), my_atom_ipi())       # (i pi)^2
def my_atom_pi4():   return mul(my_atom_pi2(), my_atom_pi2())       # (i pi)^4
