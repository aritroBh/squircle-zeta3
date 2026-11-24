# g9_r2_altmzv.py -- own exact re-derivation of the 13 values H_w(-1), weight <= 5, in the 20 constants of Section 32 (edge r = 2).
# Independent of the notes: no import of, and no data from, the notes' code or pickles. The only input read from this folder is the
# text of g9_hex5_r2fold.py (parsed with ast, not executed): the 13 words and the 13 values that the section uses now.
#
# Alternating multiple zeta values (level 2):  zeta(s; eps) = sum_{n1>...>nk>=1} prod eps_i^{n_i} n_i^{-s_i},  eps_i = +-1,
# convergent iff (s1, eps1) != (1, 1).  Iterated integrals  G(a1..aN) = int_{1>t1>...>tN>0} prod dt_i/(t_i - a_i),  a_i in {0,1,-1}.
# Dictionary (derived in the section): zeta(s; eps) = (-1)^k G(w),  w = (0^{s1-1}, e1, 0^{s2-1}, e2, ..., 0^{sk-1}, ek),  e_j = eps_1...eps_j.
# Convergent words: first letter != 1, last letter != 0 (they correspond one to one to convergent (s; eps)).
#
# Relation families (each is a proved theorem; all are applied with exact rational arithmetic, never a numerical fit):
#  SH   shuffle:   G(u)G(v) = G(u sh v) for convergent words u, v. Iterated integrals over the same path (Chen 1977; Goncharov 1995, 2001).
#  ST   stuffle:   zeta(u)zeta(v) = sum over the quasi-shuffles of the nested sums, with (s,eps)+(s',eps') = (s+s', eps eps')
#                  (Hoffman 2000; Borwein-Bradley-Broadhurst-Lisonek 2001). Series convergent (absolutely if s1 >= 2; for the symbol (1,-1), conditionally, by the limit of the finite-cutoff identity).
#                  SH - ST (transported by the dictionary) is the finite part of the coloured double shuffle relations (Racinet 2002;
#                  Ihara-Kaneko-Zagier 2006 for the regularisation theorem, which we do NOT need: only convergent words are used).
#  DUAL duality:   G(w) = (-1)^N G(tau(a_N),...,tau(a_1)),  tau(0) = (1)-(-1), tau(1) = (0)-(-1), tau(-1) = -(-1)
#                  (change of variables t -> (1-t)/(1+t) in the iterated integral; Borwein et al. 2001; Zhao 2010).
#  DIST distribution (N=2):  sum_{eta in {+-1}^k} zeta(s; eta) = 2^(k-|s|) zeta(s; 1^k) for s1 >= 2 (absolutely convergent sums,
#                  keep the even n_i; Racinet 2002, Deligne-Goncharov 2005). Used as a consistency check: it must not change the ranks.
# Not used: regularised (divergent-word) relations and the regularisation comparison; they are not needed (see the ranks printed).
import ast, os, sys, time, itertools
from fractions import Fraction as Fr
from collections import defaultdict
import sympy as sp, mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
T0 = time.time()

# ------------------------------------------------------------------ words, symbols, dictionary
def sym_to_word(s, eps):
    w, e = [], 1
    for si, ei in zip(s, eps):
        e *= ei
        w += [0]*(si-1) + [e]
    return tuple(w)

def word_to_sym(w):
    s, eps, run, e = [], [], 0, 1
    for a in w:
        if a == 0: run += 1
        else:
            s.append(run+1); eps.append(a*e); e = a; run = 0   # a = e_j, eps_j = e_j / e_{j-1} = e_j * e_{j-1}
    assert run == 0
    return tuple(s), tuple(eps)

def conv_words(n):
    out = []
    for w in itertools.product((0, 1, -1), repeat=n):
        if w[0] != 1 and w[-1] != 0: out.append(w)
    return out

def depth(w): return sum(1 for a in w if a != 0)

def vadd(d, w, c):
    d[w] = d.get(w, 0) + c
    if d[w] == 0: del d[w]

def shuffle(u, v):
    res = defaultdict(int)
    def rec(u, v, pre):
        if not u: res[pre+v] += 1; return
        if not v: res[pre+u] += 1; return
        rec(u[1:], v, pre+(u[0],)); rec(u, v[1:], pre+(v[0],))
    rec(u, v, ())
    return res

def stuffle(u, v):               # u, v: tuples of (s,eps) pairs; returns dict symbol -> multiplicity
    res = defaultdict(int)
    def rec(u, v, pre):
        if not u: res[pre+v] += 1; return
        if not v: res[pre+u] += 1; return
        rec(u[1:], v, pre+(u[0],)); rec(u, v[1:], pre+(v[0],))
        rec(u[1:], v[1:], pre+((u[0][0]+v[0][0], u[0][1]*v[0][1]),))
    rec(u, v, ())
    return res

def zeta_vec(sym):               # zeta(symbol) as a vector over words: (-1)^k G(w)
    s = tuple(p[0] for p in sym); eps = tuple(p[1] for p in sym)
    return {sym_to_word(s, eps): (-1)**len(s)}

TAU = {0: {(1,): 1, (-1,): -1}, 1: {(0,): 1, (-1,): -1}, -1: {(-1,): -1}}
def dual(w):                     # (-1)^N G(tau(a_N) ... tau(a_1)) as a vector over words
    res = {(): 1}
    for a in reversed(w):
        new = defaultdict(int)
        for pre, c in res.items():
            for t, ct in TAU[a].items(): new[pre+t] += c*ct
        res = new
    return {k: Fr((-1)**len(w)*c) for k, c in res.items() if c != 0}

# ------------------------------------------------------------------ relations by weight
def relations(n, fams):
    C = conv_words(n); idx = {w: i for i, w in enumerate(C)}
    rows = []
    def push(vec):
        row = [Fr(0)]*len(C)
        for w, c in vec.items(): row[idx[w]] += c           # KeyError would mean a non-convergent word: a bug
        if any(row): rows.append(row)
    if 'SH' in fams or 'ST' in fams:
        for a in range(1, n//2 + 1):
            b = n - a
            for u in conv_words(a):
                for v in conv_words(b):
                    sh = {w: Fr(c) for w, c in shuffle(u, v).items()}
                    su, sv = word_to_sym(u), word_to_sym(v)
                    pu = tuple(zip(*su)); pv = tuple(zip(*sv))
                    st = defaultdict(Fr)
                    for sym, m in stuffle(pu, pv).items():
                        for w, c in zeta_vec(sym).items(): st[w] += m*c
                    sgn = (-1)**(len(su[0]) + len(sv[0]))
                    # G(u)G(v) = (-1)^(ku+kv) zeta(u) zeta(v)
                    rel = defaultdict(Fr)
                    for w, c in sh.items(): rel[w] += c
                    for w, c in st.items(): rel[w] -= sgn*c
                    push({w: c for w, c in rel.items() if c != 0})
    if 'DUAL' in fams:
        for w in C:
            rel = defaultdict(Fr); rel[w] += 1
            for x, c in dual(w).items(): rel[x] -= c
            push({x: c for x, c in rel.items() if c != 0})
    if 'DIST' in fams:
        for k in range(1, n+1):
            for comp in itertools.product(range(1, n+1), repeat=k):
                if sum(comp) != n or comp[0] < 2: continue
                rel = defaultdict(Fr)
                for eta in itertools.product((1, -1), repeat=k):
                    for w, c in zeta_vec(tuple(zip(comp, eta))).items(): rel[w] += c
                for w, c in zeta_vec(tuple((s, 1) for s in comp)).items(): rel[w] -= Fr(2)**(k-n)*c
                push({w: c for w, c in rel.items() if c != 0})
    return C, idx, rows

# ------------------------------------------------------------------ exact linear algebra over Q (Fractions)
def rref(rows, ncols):
    rows = [r[:] for r in rows]; piv = []; r0 = 0
    for c in range(ncols):
        p = next((i for i in range(r0, len(rows)) if rows[i][c] != 0), None)
        if p is None: continue
        rows[r0], rows[p] = rows[p], rows[r0]
        inv = 1/rows[r0][c]; rows[r0] = [x*inv for x in rows[r0]]
        for i in range(len(rows)):
            if i != r0 and rows[i][c] != 0:
                f = rows[i][c]; ri = rows[r0]
                rows[i] = [x - f*y for x, y in zip(rows[i], ri)]
        piv.append(c); r0 += 1
        if r0 == len(rows): break
    return rows[:r0], piv

def normal_form(vec, R, piv):
    vec = vec[:]
    for row, c in zip(R, piv):
        if vec[c] != 0:
            f = vec[c]; vec = [x - f*y for x, y in zip(vec, row)]
    return vec

def solve_coeffs(target, basis):          # target = sum c_i basis_i exactly (vectors over Q); returns c or None
    m = len(basis); n = len(target)
    # unknowns c_1..c_m: columns of the system are the basis vectors; augmented column is the target
    A = [[basis[j][i] for j in range(m)] + [target[i]] for i in range(n)]
    Rr, piv = rref(A, m + 1)
    if m in piv: return None              # inconsistent
    c = [Fr(0)]*m
    for row, p in zip(Rr, piv): c[p] = row[m]
    if len(piv) < m: raise ValueError('named vectors dependent modulo the relations')
    return c

# ------------------------------------------------------------------ the 20 named constants as vectors over convergent words
def vmul(a, b):                            # product of two G-vectors by the shuffle product
    res = defaultdict(Fr)
    for u, cu in a.items():
        for v, cv in b.items():
            for w, m in shuffle(u, v).items(): res[w] += cu*cv*m
    return {w: c for w, c in res.items() if c != 0}
def vscale(a, s): return {w: c*s for w, c in a.items()}
def vsum(*vs):
    res = defaultdict(Fr)
    for v in vs:
        for w, c in v.items(): res[w] += c
    return {w: c for w, c in res.items() if c != 0}
def G1(*w): return {tuple(w): Fr(1)}

l2v = G1(-1)                               # zeta(1bar) = -log 2 = -G(-1)  =>  log 2 = G(-1)
z2v = vscale(G1(0, 1), -1)                 # zeta(2) = -G(0,1)
z3v = vscale(G1(0, 0, 1), -1)
z4v = vscale(G1(0, 0, 0, 1), -1)
z5v = vscale(G1(0, 0, 0, 0, 1), -1)
pi2v = vscale(z2v, 6); pi4v = vscale(z4v, 90)
def Lih(n):                                # Li_n(1/2): t = s/(1+s):  dt/t -> ds/s - ds/(1+s),  dt/(1-t) -> ds/(1+s)
    res = defaultdict(Fr)
    for pick in itertools.product((0, -1), repeat=n-1):
        res[tuple(pick) + (-1,)] += (-1)**sum(1 for a in pick if a == -1)
    return dict(res)
Li4v, Li5v = Lih(4), Lih(5)
def P(*fs):
    r = fs[0]
    for f in fs[1:]: r = vmul(r, f)
    return r
def pw(f, k): return P(*([f]*k))
sym = dict(zip(['l2', 'pi', 'z3', 'z5', 'Li4h', 'Li5h'], sp.symbols('l2 pi z3 z5 Li4h Li5h')))
l2, pi, z3, z5, Li4h, Li5h = [sym[k] for k in ['l2', 'pi', 'z3', 'z5', 'Li4h', 'Li5h']]
# (name, sympy value, G-vector, weight)  -- the order and the names of the section (g9_hex5_r2fold.py, `names`)
NAMED = [
 ('l2', l2, l2v, 1),
 ('pi2', pi**2, pi2v, 2), ('l2^2', l2**2, pw(l2v, 2), 2),
 ('z3', z3, z3v, 3), ('pi2 l2', pi**2*l2, P(pi2v, l2v), 3), ('l2^3', l2**3, pw(l2v, 3), 3),
 ('pi4', pi**4, pi4v, 4), ('z3 l2', z3*l2, P(z3v, l2v), 4), ('pi2 l2^2', pi**2*l2**2, P(pi2v, pw(l2v, 2)), 4),
 ('l2^4', l2**4, pw(l2v, 4), 4), ('Li4h', Li4h, Li4v, 4),
 ('z5', z5, z5v, 5), ('z2z3', pi**2*z3/6, P(z2v, z3v), 5), ('pi4 l2', pi**4*l2, P(pi4v, l2v), 5),
 ('z3 l2^2', z3*l2**2, P(z3v, pw(l2v, 2)), 5), ('pi2 l2^3', pi**2*l2**3, P(pi2v, pw(l2v, 3)), 5),
 ('l2^5', l2**5, pw(l2v, 5), 5), ('Li5h', Li5h, Li5v, 5), ('Li4h l2', Li4h*l2, P(Li4v, l2v), 5),
]

# ------------------------------------------------------------------ the 13 words and the values used by the section (read, not executed)
src = open(os.path.join(HERE, 'g9_hex5_r2fold.py')).read()
tree = ast.parse(src)
def grab(name):
    for nd in tree.body:
        if isinstance(nd, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name for t in nd.targets):
            return nd.value
    raise KeyError(name)
ns = {'l2': l2, 'pi': pi, 'z3': z3, 'z5': z5, 'Li4h': Li4h, 'Li5h': Li5h, 'R': sp.Rational, 'Z23': lambda c: c*pi**2*z3/6}
atomval = eval(compile(ast.Expression(grab('atomval')), 'atomval', 'eval'), ns)
NAMEWORD = ast.literal_eval(grab('NAMEWORD'))
assert len(atomval) == 13 and set(atomval) == set(NAMEWORD)
def hword_to_G(hw):                         # H_w(-1) = (-1)^{#1} G(w[1 -> -1]);  derived from the nested series (see below)
    return tuple(-1 if ch == '1' else 0 for ch in hw), (-1)**hw.count('1')
# H_w(x), w = 0^{s1-1}1 0^{s2-1}1 ... = sum_{n1>...>nk} x^{n1}/(n1^{s1}...nk^{sk}) = zeta(s; (-1,1,...,1)) at x=-1, whose e_j are all -1.
def hword_symbol(hw):
    s, run = [], 0
    for ch in hw:
        if ch == '0': run += 1
        else: s.append(run+1); run = 0
    return tuple(s)
for nm, hw in NAMEWORD.items():                # the dictionary reproduces the word H_w -> G(w[1->-1])
    s = hword_symbol(hw); eps = (-1,) + (1,)*(len(s)-1)
    assert sym_to_word(s, eps) == hword_to_G(hw)[0] and (-1)**len(s) == hword_to_G(hw)[1]

# ------------------------------------------------------------------ run: ranks per weight and family
print('Own exact reduction of alternating MZVs, weight <= 5 (relations: SH, ST, DUAL, DIST; exact rational arithmetic)')
print('weight | convergent words | rank SH+ST | rank DUAL | rank SH+ST+DUAL | rank +DIST | codim | Deligne bound (Fibonacci)')
data = {}; FIB = {1: 1, 2: 2, 3: 3, 4: 5, 5: 8}
for n in range(1, 6):
    ranks = {}
    for tag, fams in [('DS', ('SH', 'ST')), ('DUAL', ('DUAL',)), ('DS+DUAL', ('SH', 'ST', 'DUAL')), ('ALL', ('SH', 'ST', 'DUAL', 'DIST'))]:
        C, idx, rows = relations(n, fams)
        R, piv = rref(rows, len(C)) if rows else ([], [])
        ranks[tag] = (len(rows), len(R))
        if tag == 'ALL': dataall = (C, idx, R, piv)
    print('  %d    |   %3d            |   %3d      |   %3d     |   %3d           |   %3d      |  %d   |  %d' % (
        n, len(C), ranks['DS'][1], ranks['DUAL'][1], ranks['DS+DUAL'][1], ranks['ALL'][1], len(C)-ranks['ALL'][1], FIB[n]))
    data[n] = dataall                            # reduce modulo all four families
    assert len(C) - ranks['ALL'][1] == FIB[n], 'codimension is not the Deligne bound'
print('relations counted before elimination (SH+ST, DUAL, DIST):')
for n in range(1, 6):
    print('  weight %d:' % n, {f: len(relations(n, (f,) if f != 'SH+ST' else ('SH', 'ST'))[2]) for f in ['SH+ST', 'DUAL', 'DIST']})
print('(%.0f s)' % (time.time()-T0))

# ------------------------------------------------------------------ reduce the named constants and the 13 targets modulo the relations
def to_row(vec, C, idx):
    row = [Fr(0)]*len(C)
    for w, c in vec.items(): row[idx[w]] += c
    return row
def coeffs_in_named(vec, n):
    C, idx, R, piv = data[n]
    named = [(nm, e, v) for nm, e, v, wt in NAMED if wt == n]
    NF = [normal_form(to_row(v, C, idx), R, piv) for _, _, v in named]
    t = normal_form(to_row(vec, C, idx), R, piv)
    c = solve_coeffs(t, NF)
    return named, c, NF
print('\nNamed constants modulo the relations (independent in each weight, hence a basis of the quotient of dimension = codim):')
for n in range(1, 6):
    named, _, NF = coeffs_in_named({}, n)
    A, piv = rref(NF, len(NF[0]))
    print('  weight %d: %d named constants, rank of their images in the quotient = %d, quotient dimension = %d' % (n, len(named), len(piv), FIB[n]))
    assert len(piv) == len(named) == FIB[n]

derived = {}
print('\nDerived values H_w(-1) (every one is an exact solution of  T - sum c_i N_i in span(relations)):')
for nm, hw in NAMEWORD.items():
    gw, sg = hword_to_G(hw)
    n = len(hw)
    named, c, NF = coeffs_in_named({gw: Fr(sg)}, n)
    if c is None:
        print('  %-7s word %-6s NOT DETERMINED by the relations generated' % (nm, hw)); sys.exit(1)
    # certificate check: remainder reduces to zero
    C, idx, R, piv = data[n]
    resid = normal_form(to_row({gw: Fr(sg)}, C, idx), R, piv)
    for ci, f in zip(c, NF): resid = [x - ci*y for x, y in zip(resid, f)]
    assert not any(resid)
    val = sum(sp.Rational(ci.numerator, ci.denominator)*e for ci, (_, e, _) in zip(c, named))
    derived[nm] = sp.expand(val)
    print('  %-7s = %s' % (nm, derived[nm]))

# ------------------------------------------------------------------ exact comparison with the values used by the section
print('\nExact comparison (derived - value used in g9_hex5_r2fold.py), expanded over the polynomial ring Q[l2, pi, z3, z5, Li4h, Li5h]:')
allzero = True
for nm in NAMEWORD:
    d = sp.expand(derived[nm] - sp.expand(atomval[nm]))
    print('  %-7s difference = %s' % (nm, d)); allzero &= (d == 0)
print('all 13 differences are exactly 0:', allzero)
assert allzero

# ------------------------------------------------------------------ independent numerical check, mpmath, 80 digits
mp.mp.dps = 90
def nested_terms(ss, N):                    # b_n = n^{-s1} sum_{n>n2>...>nk} prod n_i^{-s_i}, n = 1..N
    arr = [mp.mpf(1)]*(N+2)
    for s in reversed(ss[1:]):
        new = [mp.mpf(0)]*(N+2); acc = mp.mpf(0)
        for j in range(1, N+1):
            acc += arr[j]/mp.mpf(j)**s
            new[j+1] = acc
        arr = new
    return [None] + [arr[n]/mp.mpf(n)**ss[0] for n in range(1, N+1)]
def cvz(a, N):                              # Cohen-Villegas-Zagier, Algorithm 1: sum_{k>=0} (-1)^k a_k from a_0..a_{N-1}
    d = (3 + mp.sqrt(8))**N; d = (d + 1/d)/2; b = mp.mpf(-1); c = -d; s = mp.mpf(0)
    for k in range(N):
        c = b - c; s += c*a[k]
        b = (k+N)*(k-N)*b/((k+mp.mpf(1)/2)*(k+1))
    return s/d
def Hm1(hw, N):                             # H_w(-1) = sum_{n>=1} (-1)^n b_n
    b = nested_terms(hword_symbol(hw), N+2)
    return -cvz([b[k+1] for k in range(N)], N)
cv = {l2: mp.log(2), pi: mp.pi, z3: mp.zeta(3), z5: mp.zeta(5), Li4h: mp.polylog(4, mp.mpf(1)/2), Li5h: mp.polylog(5, mp.mpf(1)/2)}
def ev(e): return sp.lambdify(list(cv), e, 'mpmath')(*cv.values())
print('\nNumerical check (mpmath, 90 digits; direct nested series summed with Cohen-Villegas-Zagier acceleration, 130 and 150 terms):')
worst = mp.mpf(0); worstN = mp.mpf(0)
for nm, hw in NAMEWORD.items():
    a = Hm1(hw, 130); b = Hm1(hw, 150); d = abs(a - ev(derived[nm]))
    worst = max(worst, d); worstN = max(worstN, abs(a-b))
    print('  %-7s series = %s   |series - derived| = %s   (130 vs 150 terms: %s)' % (nm, mp.nstr(b, 55), mp.nstr(d, 3), mp.nstr(abs(a-b), 3)))
print('worst |series - derived| over the 13 values = %s ;  worst change between 130 and 150 terms = %s' % (mp.nstr(worst, 3), mp.nstr(worstN, 3)))
assert worst < mp.mpf(10)**-60
print('total time %.0f s' % (time.time()-T0))
