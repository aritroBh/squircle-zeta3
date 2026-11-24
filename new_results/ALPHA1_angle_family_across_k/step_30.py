"""step_30.py -- exact HPL(c) form of J(r) = slope*(2c^2/(1-c^2)) via a04's GPL->HPL machinery.
J = int_0^1 [ell^2 (f_P - f_Q2) + ell emm (f_Q - 2 f_P)] dx, letters P=-1/c, Q=-(1-c)/(2c), Q2=-(1-c)/c.
Numeric check vs quadrature J(r) at several c."""
import sys
sys.path.insert(0, '.')
from fractions import Fraction as Fr
from collections import defaultdict
import mpmath as mp
mp.mp.dps = 50
from hpl_tools import shuffle_words, clean

SAMEPOS = {('E', '1'), ('1', 'E'), ('B', '0'), ('0', 'B'), ('E', 'E'), ('B', 'B')}
def canon(x):
    return {'E': '1', 'B': '0'}.get(x, x)
DL = {}
def setd(x, y, d):
    DL[(x, y)] = d; DL[(y, x)] = d
setd('1', '0', {})
setd('1', 'P', {-1: 1, 0: -1})
setd('1', 'Q', {-1: 1, 0: -1})
setd('1', 'Q2', {0: -1})
setd('0', 'P', {0: -1})
setd('0', 'Q', {1: -1, 0: -1})
setd('0', 'Q2', {1: -1, 0: -1})
setd('P', 'Q', {-1: 1, 0: -1})
setd('P', 'Q2', {})
setd('Q', 'Q2', {1: -1, 0: -1})
def D(x, y):
    x, y = canon(x), canon(y)
    if x == y: return {}
    return DL[(x, y)]

ZETA = {}
def mzv_const(word):
    if word in ZETA: return ZETA[word]
    if len(word) == 0:
        res = {(): Fr(1)}
    elif all(a == '0' for a in word):
        res = {}
    elif word[0] == '1':
        raise ValueError('divergent constant word: ' + str(word))
    elif word[-1] == '0':
        m = 0
        while word[-1 - m] == '0': m += 1
        base = word[:-1]
        acc = defaultdict(Fr)
        for w, cnt in shuffle_words(base, ('0',)).items():
            if w == word: continue
            for mono, v in mzv_const(w).items():
                acc[mono] -= cnt*v
        res = {k: v/m for k, v in acc.items() if v}
    else:
        table = {('0', '1'): {('z2',): Fr(-1)},
                 ('0', '0', '1'): {('z3',): Fr(-1)},
                 ('0', '1', '1'): {('z3',): Fr(1)},
                 ('0', '0', '0', '1'): {('z4',): Fr(-1)},
                 ('0', '1', '0', '1'): {('z2','z2'): Fr(1,4)},
                 ('0', '0', '1', '1'): {('z4',): Fr(1,2)},
                 ('0', '1', '1', '1'): {('z4',): Fr(-1)}}
        # z4 = zeta(4) = pi^4/90; 0,1,0,1 = Li4(1/2)-... use z2^2 coeff 1/4 (z(4) = 4/90 pi^4?) -- fill below
        if word not in table:
            raise ValueError('MZV value not tabulated: ' + str(word))
        res = table[word]
    ZETA[word] = res
    return res

MEMO = {}
def conv(word):
    if word in MEMO: return MEMO[word]
    if all(a in ('0', '1') for a in word):
        res = {((), mono): v for mono, v in mzv_const(word).items()}
        MEMO[word] = res
        return res
    n = len(word)
    pts = ('E',) + tuple(word) + ('B',)
    acc = defaultdict(Fr)
    for i in range(1, n + 1):
        form = defaultdict(Fr)
        for a, v in D(pts[i - 1], pts[i]).items(): form[a] += v
        for a, v in D(pts[i + 1], pts[i]).items(): form[a] -= v
        form = {a: v for a, v in form.items() if v}
        if not form: continue
        sub = conv(word[:i - 1] + word[i:])
        for a, fv in form.items():
            for (w, mono), sv in sub.items():
                acc[((a,) + w, mono)] += fv*sv
    res = {k: v for k, v in acc.items() if v}
    MEMO[word] = res
    return res

def gshuffle(A, B):
    out = defaultdict(Fr)
    for (u, ku), a in A.items():
        for (v, kv), b in B.items():
            for w, cnt in shuffle_words(u, v).items():
                out[(w, ku + kv)] += a*b*cnt
    return {k: v for k, v in out.items() if v}

ell = {(('0',), 0): Fr(1), (('P',), 0): Fr(1), (('Q',), 0): Fr(-1), ((), 1): Fr(1)}
emm = {(('1',), 0): Fr(1), (('Q2',), 0): Fr(1), (('Q',), 0): Fr(-1)}
ell2 = gshuffle(ell, ell); ellm = gshuffle(ell, emm)

def outer(expr, letter, coef):
    return {((letter,) + w, k): coef*v for (w, k), v in expr.items()}

Jexpr = defaultdict(Fr)
for k, v in outer(ell2, 'P', 1).items(): Jexpr[k] += v
for k, v in outer(ell2, 'Q2', -1).items(): Jexpr[k] += v
for k, v in outer(ellm, 'Q', 1).items(): Jexpr[k] += v
for k, v in outer(ellm, 'P', -2).items(): Jexpr[k] += v

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
Jc = to_hpl(Jexpr)
print("J(c): %d GPL terms -> %d HPL terms" % (len(Jexpr), len(Jc)))
monos = defaultdict(int)
for (w, mono) in Jc: monos[mono] += 1
print("monomials:", dict(monos))

# numeric evaluation: H_w(c) via series; constants z2, z3, z4
import hpl_tools as HT
def hpl_num(word, c, N=60):
    # Taylor series in c of H_word(c)
    ser = HT.hpl_taylor(word, N)
    return mp.polyval([mp.mpf(int(x.numerator))/int(x.denominator) for x in ser], c)
def mono_val(mono, c):
    val = mp.mpf(1)
    for m in mono:
        if m == 'z2': val *= mp.zeta(2)
        elif m == 'z3': val *= mp.zeta(3)
        elif m == 'z4': val *= mp.zeta(4)
    return val
def J_hpl(c):
    tot = mp.mpf(0)
    for (w, mono), v in Jc.items():
        tot += mp.mpf(int(v.numerator))/int(v.denominator) * hpl_num(w, c) * mono_val(mono, c)
    return tot

def Jval(c):
    r = (1-c)/(1+c)
    def lt(w): return mp.log((w-r)*(w+1)/(2*(1-r)*w))
    def lmt(w): return mp.log((1-w)*(w+r)/(2*(1-r)*w))
    def f(w):
        return (1/(w+1) - 1/(w+r))*lt(w)**2 + (1/w - 2/(w+1))*lt(w)*lmt(w)
    return mp.quad(f, [r, mp.sqrt(r), 1])

for c in (mp.mpf(1)/2, mp.mpf(3)/5, mp.mpf(1)/3, mp.mpf(4)/5, mp.mpf(2)/7):
    Jh = J_hpl(c); Jn = Jval(c)
    print(f"c={mp.nstr(c,8)}: J_hpl={mp.nstr(Jh,25)}  J_quad={mp.nstr(Jn,25)}  diff={mp.nstr(abs(Jh-Jn),3)}")
print("STEP30 DONE")
