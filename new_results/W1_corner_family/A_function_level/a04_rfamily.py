"""a04_rfamily.py -- EXACT HPL(c) representation of the squircle pieces U, V (A3 split W4 = U + sqrt2 V).

Scale invariance.  With t(v) = (v-b)(v+a)/(4v) (the paper's t), put v = a w.  Since 4 = 2(a-b) and b/a = b^2,
   t = (w - r)(w + 1) / (2(1 - r) w),   1 - t = (1 - w)(w + r) / (2(1 - r) w),   r = b/a = b^2 = 3 - 2 sqrt2,
dv/v = dw/w, dv/(v+a) = dw/(w+1).  So U = U_r, V = V_r at r = b^2, where for 0 < r < 1
   U_r = int_r^1 l^3 (1/(6w) - 1/(3(w+1))) dw,     V_r = int_r^1 (-7 l^3/(48 w) + l^2 m/(16 w)) dw,
   l = log t, m = log(1-t).
With w = r + (1-r) x (x in [0,1]) and c := (1-r)/(1+r)  (so r = b^2  <=>  c = 1/sqrt2):
   l = G(0;x) + G(P;x) - G(Q;x) + H_1(c),   m = G(1;x) + G(Q2;x) - G(Q;x),
   P = -1/c, Q = -(1-c)/(2c), Q2 = -(1-c)/c  (positions of the GPL letters),
   dw/w = dx/(x-Q),  dw/(w+1) = dx/(x-P).
[check: log(x+p) - log(x+q) - log 2 = G(P;x) - G(Q;x) + log(p/(2q)), p/(2q) = 1/(1-c);
        log(1-x) + log(x+2q) - log(x+q) - log 2 = G(1;x) + G(Q2;x) - G(Q;x).]
Every letter difference has dlog in span{dlog c, dlog(1-c), dlog(1+c)}, and every non-constant GPL -> 0 as
c -> 0 (letters P, Q, Q2 -> -infinity).  Hence U_r, V_r are exact Q[zeta2,zeta3]-combinations of HPLs H_w(c)
with nonzero last letter.  This script derives them exactly, then checks numerically against direct
quadrature of the squircle integrals U, V at c = 1/sqrt2.
"""
import sys, json, pickle
from fractions import Fraction as Fr
from collections import defaultdict
import mpmath as mp
from hpl_tools import shuffle_words, clean

# positions: 'E' endpoint 1, 'B' base 0; letters '0','1','P','Q','Q2'
SAMEPOS = {('E', '1'), ('1', 'E'), ('B', '0'), ('0', 'B'), ('E', 'E'), ('B', 'B')}
# dlog of differences, as dict over HPL letters {0:f0, 1:f1, -1:f_-1}:  dlog c = f0, dlog(1-c) = -f1, dlog(1+c) = f_-1
def canon(x):
    return {'E': '1', 'B': '0'}.get(x, x)

DL = {}
def setd(x, y, d):
    DL[(x, y)] = d; DL[(y, x)] = d
setd('1', '0', {})                       # 1 - 0 = 1
setd('1', 'P', {-1: 1, 0: -1})           # 1 + 1/c = (1+c)/c
setd('1', 'Q', {-1: 1, 0: -1})           # (1+c)/(2c)
setd('1', 'Q2', {0: -1})                 # 1/c
setd('0', 'P', {0: -1})                  # 1/c
setd('0', 'Q', {1: -1, 0: -1})           # (1-c)/(2c)
setd('0', 'Q2', {1: -1, 0: -1})          # (1-c)/c
setd('P', 'Q', {-1: 1, 0: -1})           # -(1+c)/(2c)
setd('P', 'Q2', {})                      # -1
setd('Q', 'Q2', {1: -1, 0: -1})          # (1-c)/(2c)

def D(x, y):
    x, y = canon(x), canon(y)
    if x == y:
        return {}
    return DL[(x, y)]

# Lin: dict (word, mono) -> Fraction ; mono in {(), ('z2',), ('z3',)}
def lin_add(acc, d, s=Fr(1)):
    for k, v in d.items():
        acc[k] += s*v

ZETA = {}
def mzv_const(word):
    """G(word;1) for word in {'0','1'}^*, first letter '0' (or empty), shuffle-regularised G(0;1)=0.
    Returns dict mono -> Fraction (weight <= 3 needed)."""
    if word in ZETA:
        return ZETA[word]
    if len(word) == 0:
        res = {(): Fr(1)}
    elif all(a == '0' for a in word):
        res = {}
    elif word[0] == '1':
        raise ValueError('divergent constant word needed: ' + str(word))
    elif word[-1] == '0':
        # m * G(w 0^m) = - sum_{other insertions of the last 0} G(...)
        m = 0
        while word[-1 - m] == '0':
            m += 1
        base = word[:-1]
        acc = defaultdict(Fr)
        for w, cnt in shuffle_words(base, ('0',)).items():
            if w == word:
                continue
            for mono, v in mzv_const(w).items():
                acc[mono] -= cnt*v
        res = {k: v/m for k, v in acc.items() if v}
    else:
        table = {('0', '1'): {('z2',): Fr(-1)},
                 ('0', '0', '1'): {('z3',): Fr(-1)},
                 ('0', '1', '1'): {('z3',): Fr(1)}}
        if word not in table:
            raise ValueError('MZV value not tabulated: ' + str(word))
        res = table[word]
    ZETA[word] = res
    return res

MEMO = {}
def conv(word):
    """HPL(c) expansion of G(word; 1): dict (hplword, mono) -> Fraction."""
    if word in MEMO:
        return MEMO[word]
    if all(a in ('0', '1') for a in word):
        res = {((), mono): v for mono, v in mzv_const(word).items()}
        MEMO[word] = res
        return res
    n = len(word)
    pts = ('E',) + tuple(word) + ('B',)
    acc = defaultdict(Fr)
    for i in range(1, n + 1):
        form = defaultdict(Fr)
        for a, v in D(pts[i - 1], pts[i]).items():
            form[a] += v
        for a, v in D(pts[i + 1], pts[i]).items():
            form[a] -= v
        form = {a: v for a, v in form.items() if v}
        if not form:
            continue
        sub = conv(word[:i - 1] + word[i:])
        for a, fv in form.items():
            for (w, mono), sv in sub.items():
                acc[((a,) + w, mono)] += fv*sv
    res = {k: v for k, v in acc.items() if v}
    MEMO[word] = res
    return res

# --- build the GPL(x) expansions of l and m;  'K' marks the constant-in-x factor H_1(c)
def gshuffle(A, B):
    out = defaultdict(Fr)
    for (u, ku), a in A.items():
        for (v, kv), b in B.items():
            for w, cnt in shuffle_words(u, v).items():
                out[(w, ku + kv)] += a*b*cnt
    return {k: v for k, v in out.items() if v}

ell = {(('0',), 0): Fr(1), (('P',), 0): Fr(1), (('Q',), 0): Fr(-1), ((), 1): Fr(1)}   # (word, power of K)
emm = {(('1',), 0): Fr(1), (('Q2',), 0): Fr(1), (('Q',), 0): Fr(-1)}
ell2 = gshuffle(ell, ell); ell3 = gshuffle(ell2, ell); ell2m = gshuffle(ell2, emm)

def outer(expr, letter, coef):
    """coef * int_0^1 expr(x) dx/(x - letter)  ->  dict (gplword, Kpower) -> Fraction"""
    return {((letter,) + w, k): coef*v for (w, k), v in expr.items()}

Uexpr = defaultdict(Fr); Vexpr = defaultdict(Fr)
for k, v in outer(ell3, 'Q', Fr(1, 6)).items(): Uexpr[k] += v
for k, v in outer(ell3, 'P', Fr(-1, 3)).items(): Uexpr[k] += v
for k, v in outer(ell3, 'Q', Fr(-7, 48)).items(): Vexpr[k] += v
for k, v in outer(ell2m, 'Q', Fr(1, 16)).items(): Vexpr[k] += v

# H_1(c)^k as HPL: shuffle power of the word (1,)
def hshuffle(A, B):
    out = defaultdict(Fr)
    for (u, mu), a in A.items():
        for (v, mv), b in B.items():
            for w, cnt in shuffle_words(u, v).items():
                out[(w, tuple(sorted(mu + mv)))] += a*b*cnt
    return {k: v for k, v in out.items() if v}

H1pow = {0: {((), ()): Fr(1)}}
for k in range(1, 4):
    H1pow[k] = hshuffle(H1pow[k - 1], {((1,), ()): Fr(1)})

def to_hpl(expr):
    out = defaultdict(Fr)
    for (gw, kp), v in expr.items():
        prod = hshuffle(H1pow[kp], conv(gw))
        for key, pv in prod.items():
            out[key] += v*pv
    return {k: v for k, v in out.items() if v}

Uc = to_hpl(Uexpr); Vc = to_hpl(Vexpr)
print('U_r(c): %d GPL terms -> %d (HPL word, constant) terms' % (len(Uexpr), len(Uc)))
print('V_r(c): %d GPL terms -> %d (HPL word, constant) terms' % (len(Vexpr), len(Vc)))
for name, X in (('U', Uc), ('V', Vc)):
    monos = defaultdict(int); wts = defaultdict(int)
    for (w, mono) in X:
        monos[mono] += 1; wts[(len(w), mono)] += 1
    print(name, 'constant monomials:', dict(monos), ' (weight, mono) counts:', dict(wts))
    assert all(w[-1] != 0 for (w, mono) in X if w), 'trailing zero in HPL word'
pickle.dump({'U': Uc, 'V': Vc}, open('a04_rfamily.pkl', 'wb'))

# ---------------- numerical check at c = 1/sqrt2 against direct quadrature of U, V
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 50
x0 = 1/mp.sqrt(2)
NT = int(mp.mp.dps*3.4/0.1505) + 40
from hpl_tools import hpl_numeric_all
allw = sorted({w for (w, mono) in list(Uc) + list(Vc)}, key=len)
vals = hpl_numeric_all(allw, x0, NT)
CONST = {(): mp.mpf(1), ('z2',): mp.zeta(2), ('z3',): mp.zeta(3)}
def ev(X):
    return mp.fsum(mp.mpf(v.numerator)/v.denominator*CONST[mono]*vals[w] for (w, mono), v in X.items())
Uh, Vh = ev(Uc), ev(Vc)
r2 = mp.sqrt(2); a = 1 + r2; b = r2 - 1
lt = lambda v: mp.log((v - b)*(v + a)/(4*v))
lm = lambda v: mp.log((a - v)*(v + b)/(4*v))
Uq = mp.quad(lambda v: lt(v)**3*(1/(6*v) - 1/(3*(v + a))), [b, 1, a])
Vq = mp.quad(lambda v: -7*lt(v)**3/(48*v) + lt(v)**2*lm(v)/(16*v), [b, 1, a])
print('U_r(1/sqrt2) (HPL)  =', mp.nstr(Uh, 40))
print('U (quadrature)      =', mp.nstr(Uq, 40), '  diff', mp.nstr(Uh - Uq, 5))
print('V_r(1/sqrt2) (HPL)  =', mp.nstr(Vh, 40))
print('V (quadrature)      =', mp.nstr(Vq, 40), '  diff', mp.nstr(Vh - Vq, 5))
print('W4 = U + sqrt2 V    =', mp.nstr(Uh + r2*Vh, 40))
