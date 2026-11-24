"""Referee script 2 (independent code): exact proof of the N_4, N_5 closed forms.

Own implementation of level-2 (alternating) MZVs as iterated integrals G(w;1), w over {0,1,-1},
G(a1..an;1) = int_{1>t1>..>tn>0} prod dt_i/(t_i-a_i).
Relations (each an instance of a proved identity):
  DS  : shuffle(u,v) = stuffle(u,v) for convergent words u,v (both products of absolutely
        convergent sums / integrals);
  DUAL: t -> (1-t)/(1+t):  G(a1..an) = (-1)^n G(tau(an),...,tau(a1)),
        tau(0) = (1)-(-1), tau(1) = (0)-(-1), tau(-1) = -(-1).
Series convention: zeta(s1..sd; e1..ed) = sum_{n1>..>nd} prod eps_i^{n_i}/n_i^{s_i}
  = (-1)^d G(0^{s1-1},c1,...,0^{sd-1},cd;1), c_i = eps_1...eps_i.
Then: N_k (from r1's exact symbolic formula in lam=log2 and J_{j,i}) minus the printed closed
form is shown to lie in the Q-span of the relations, weight by weight, by exact rational
elimination. Controls: zeta(bar2,1)=zeta(3)/8, DG codimensions, and perturbed (false) targets
must NOT reduce to zero.
"""
import pickle, time, itertools, sys
from fractions import Fraction as Fr
import sympy as sp

T0 = time.time()
def log(*a): print('[%6.1fs]' % (time.time() - T0), *a, flush=True)

LET = (0, 1, -1)
def conv(w): return len(w) > 0 and w[0] != 1 and w[-1] != 0
def conv_words(n): return [w for w in itertools.product(LET, repeat=n) if conv(w)]

def shuffle(u, v):
    if not u: return {v: 1}
    if not v: return {u: 1}
    out = {}
    for w, c in shuffle(u[1:], v).items():
        k = (u[0],) + w; out[k] = out.get(k, 0) + c
    for w, c in shuffle(u, v[1:]).items():
        k = (v[0],) + w; out[k] = out.get(k, 0) + c
    return out

def word_to_series(w):
    S, E, z, prev = [], [], 0, 1
    for a in w:
        if a == 0: z += 1
        else:
            S.append(z + 1); E.append(a * prev); prev = a; z = 0
    assert z == 0
    return tuple(zip(S, E))
def series_to_word(ser):
    w, c = [], 1
    for s, e in ser:
        c *= e; w += [0] * (s - 1) + [c]
    return tuple(w)

def stuffle(A, B):
    if not A: return {B: 1}
    if not B: return {A: 1}
    out = {}
    def add(k, c): out[k] = out.get(k, 0) + c
    for C, c in stuffle(A[1:], B).items(): add((A[0],) + C, c)
    for C, c in stuffle(A, B[1:]).items(): add((B[0],) + C, c)
    for C, c in stuffle(A[1:], B[1:]).items():
        add(((A[0][0] + B[0][0], A[0][1] * B[0][1]),) + C, c)
    return out

def add_to(d, k, c):
    v = d.get(k, 0) + c
    if v == 0: d.pop(k, None)
    else: d[k] = v

def ds_relation(u, v):
    r = {}
    for w, c in shuffle(u, v).items(): add_to(r, w, Fr(c))
    A, B = word_to_series(u), word_to_series(v)
    sgn = (-1) ** (len(A) + len(B))
    for C, c in stuffle(A, B).items():
        add_to(r, series_to_word(C), -Fr(sgn * (-1) ** len(C) * c))
    return r

TAU = {0: [(1, 1), (-1, -1)], 1: [(0, 1), (-1, -1)], -1: [(-1, -1)]}
def dual_relation(w):
    r = {w: Fr(1)}
    n = len(w)
    for combo in itertools.product(*[TAU[a] for a in reversed(w)]):
        word = tuple(x for x, _ in combo); c = (-1) ** n
        for _, cc in combo: c *= cc
        add_to(r, word, -Fr(c))
    return r

class Space:
    """incremental exact row echelon form over Q"""
    def __init__(self): self.piv = {}  # pivot word -> row (dict) normalised to 1 at pivot
    def order(self, w): return w
    def reduce(self, v):
        v = dict(v)
        while True:
            ks = [k for k in v if k in self.piv]
            if not ks: return v
            k = max(ks)
            c = v[k]
            for kk, cc in self.piv[k].items(): add_to(v, kk, -c * cc)
    def add(self, v):
        v = self.reduce(v)
        if not v: return False
        k = max(v); c = v[k]
        self.piv[k] = {kk: cc / c for kk, cc in v.items()}
        return True

SP = {}
for n in range(1, 6):
    sp_ = Space(); nrel = 0
    for m in range(1, n // 2 + 1):
        for u in conv_words(m):
            for v in conv_words(n - m):
                sp_.add(ds_relation(u, v)); nrel += 1
    for w in conv_words(n):
        sp_.add(dual_relation(w)); nrel += 1
    SP[n] = sp_
    log('weight %d: %d conv words, %d relations, rank %d, codim %d' % (
        n, len(conv_words(n)), nrel, len(sp_.piv), len(conv_words(n)) - len(sp_.piv)))

def in_span(vec):
    """vec: dict word->Fr, possibly mixed weight. Returns dict weight->residual."""
    byw = {}
    for w, c in vec.items(): byw.setdefault(len(w), {})[w] = c
    res = {}
    for n, v in byw.items():
        if n > 5: res[n] = v; continue
        r = SP[n].reduce(v)
        if r: res[n] = r
    return res

# word-vector algebra
def vmul(a, b):
    out = {}
    for u, cu in a.items():
        for v, cv in b.items():
            for w, c in shuffle(u, v).items(): add_to(out, w, cu * cv * c)
    return out
def vadd(a, b, cb=1):
    out = dict(a)
    for w, c in b.items(): add_to(out, w, c * cb)
    return out
ONE = {(): Fr(1)}
def vpow(a, m):
    out = ONE
    for _ in range(m): out = vmul(out, a)
    return out
LAM = {(-1,): Fr(1)}                 # log 2
PI2 = {(0, 1): Fr(-6)}               # pi^2 = 6 zeta(2) = -6 G(0,1)
Z3 = {(0, 0, 1): Fr(-1)}
Z5 = {(0, 0, 0, 0, 1): Fr(-1)}
import math
def Jw(j, i):  # J_{j,i} = j! i! G(0^{i+1}, (-1)^j)
    return {(0,) * (i + 1) + (-1,) * j: Fr(math.factorial(j) * math.factorial(i))}

# controls
r = in_span(vadd({(0, -1, -1): Fr(1)}, Z3, -Fr(1, 8)))
log('control zeta(bar2,1) = G(0,-1,-1) = zeta(3)/8 :', 'PROVED' if not r else 'FAIL %s' % r)
r = in_span(vadd(Jw(2, 0), Z3, -Fr(1, 4)))
log('control J_{2,0} = zeta(3)/4 :', 'PROVED' if not r else 'FAIL')
# J_{1,k} = k! eta(k+2): eta(n) = -zeta(bar n) = G(0^{n-1},-1)
for k in range(4):
    eta = {(0,) * (k + 1) + (-1,): Fr(1)}  # = eta(k+2)? check: zeta(bar n) = -G(0^{n-1},-1), eta = -zeta(bar n)
    r = in_span(vadd(Jw(1, k), eta, -Fr(math.factorial(k))))
    log('J_{1,%d} = %d! eta(%d) (word identity) :' % (k, k, k + 2), 'OK' if not r else 'FAIL')

# ---- N_k from r1 ----
D = pickle.load(open('r1_Nsym.pkl', 'rb'))
lam = sp.Symbol('lam')
def sym_to_vec(expr):
    out = {}
    for term in sp.Add.make_args(sp.expand(expr)):
        c, rest = term.as_coeff_Mul()
        c = sp.Rational(c); v = {(): Fr(int(c.p), int(c.q))}
        for fac in sp.Mul.make_args(rest):
            b, ex = fac.as_base_exp()
            if b == 1: continue
            if b == lam: v = vmul(v, vpow(LAM, int(ex)))
            else:
                _, j, i = str(b).split('_'); assert ex == 1
                v = vmul(v, Jw(int(j), int(i)))
        out = vadd(out, v)
    return out

P2s, Z3s, Z5s = sp.symbols('P2 Z3 Z5')
CF = {
 0: -lam, 1: lam**2 / 2, 2: P2s / 24 - lam**3 / 6,
 3: lam**4 / 24 + P2s / 24 * (1 - lam) + Z3s / 4,
 4: -lam**5 / 120 + P2s / 24 * (1 - lam + lam**2 / 2) + Z3s / 4 * (1 - lam) + P2s**2 / 2880,
 5: lam**6 / 720 + P2s / 24 * (1 - lam + lam**2 / 2 - lam**3 / 6) + Z3s / 4 * (1 - lam + lam**2 / 2)
    + P2s**2 / 2880 * (1 - lam) - 7 * P2s**2 / 1920 - Z3s * P2s / 48 + Z5s / 8,
}
def cf_to_vec(expr):
    out = {}
    for term in sp.Add.make_args(sp.expand(expr)):
        c, rest = term.as_coeff_Mul()
        c = sp.Rational(c); v = {(): Fr(int(c.p), int(c.q))}
        for fac in sp.Mul.make_args(rest):
            b, ex = fac.as_base_exp(); ex = int(ex)
            base = {lam: LAM, P2s: PI2, Z3s: Z3, Z5s: Z5}[b]
            v = vmul(v, vpow(base, ex))
        out = vadd(out, v)
    return out

for k in range(6):
    Nv = sym_to_vec(D['N'][k])
    diff = vadd(Nv, cf_to_vec(CF[k]), -1)
    r = in_span(diff)
    wts = sorted(set(len(w) for w in Nv))
    log('N_%d: weights present %s; N_k - closed form in span(DS+DUAL): %s' % (
        k, wts, 'YES (exact)' if not r else 'NO, residual at weights %s' % sorted(r)))
# negative controls: perturb closed forms
for k, pert in [(4, P2s**2 / 10**6), (5, Z5s / 10**6), (5, Z3s * P2s / 10**6), (3, Z3s / 1000)]:
    diff = vadd(sym_to_vec(D['N'][k]), cf_to_vec(CF[k] + pert), -1)
    r = in_span(diff)
    log('negative control N_%d with perturbation %s: %s' % (k, pert, 'correctly NOT in span' if r else 'WRONGLY in span'))
log('done')
