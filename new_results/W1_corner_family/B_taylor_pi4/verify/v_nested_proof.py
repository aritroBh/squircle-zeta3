"""v_nested_proof.py -- verifier's independent all-n proof of G1's B3 (and B2): exact nested-sum algebra.

Written from scratch (does not import the solver's zsum.py / hpl.py). Coefficients are exact rational
functions of n built on python-flint fmpq_poly (own class RF below), not sympy.

Objects.  Z_w(M) = sum_{M >= i1 > ... > ik >= 1} prod_j x_{w_j}(i_j),  x_a(i) = sgn(a)^i / i^|a|.
An expression is {word: RF}, meaning sum RF(n) Z_word(2n).

Facts used (each an identity valid for all M >= 1, re-checked exactly below at many n):
  (Q)  Z_{a,u} Z_{b,v} = Z_{a,(u*Z_{b,v})} + Z_{b,(Z_{a,u}*v)} + Z_{a o b,(u*v)},  a o b = sgn(a)sgn(b)(|a|+|b|).
  (P)  Z_{a,u}(M) = Z_{a,u}(M-1) + x_a(M) Z_u(M-1)   (peeling; used to move 2n+d -> 2n).
  (H)  HPL Taylor coefficients: H_w(x) = sum_N c_w(N) x^N with
          c_(1)(N) = 1/N,  c_(-1)(N) = -(-1)^N/N,  c_(0,w)(N) = c_w(N)/N,
          c_(a,w)(N) = (1/N) sum_{i<N} a^(N-1-i) c_w(i)   (a = +-1; from f_a = 1/(1 - a t)).
       If c_w(i) = k sigma^i i^-m Z_v(i-1), then c_(a,w)(N) = a k a^N N^-1 Z_{(a sigma) m, v}(N-1).
  (h)  H^(r)_{2n-1} = Z_r(2n) - (2n)^-r,   H^(r)_{n-1} = 2^(r-1)(Z_r(2n) + Z_{-r}(2n)) - n^-r.
Claim checked: for every part (rational, zeta2, zeta3) of B3 and B2, as identities of {word: RF}:
   [c^(2n)] B - mu_{k,n} == 0  and  [c^(2n+1)] B == 0   for all n >= 1,
plus the finitely many remaining coefficients (c^-2 .. c^1) in exact arithmetic.
mu_{k,n} is taken from the verifier's own derivation (v_mu_formula.pkl, from v_mu3_derive.py).
"""
import json, re, sys, pickle, random
from fractions import Fraction as Fr
from functools import lru_cache
import flint
import sympy as sp

sys.setrecursionlimit(20000)
X = flint.fmpq_poly([0, 1])          # the variable n
G1 = 'research/10_frontier/G1_zeta3_any_angle/'


class RF:
    """exact rational function num/den in Q(n), den monic, gcd(num, den) = 1."""
    __slots__ = ('p', 'q')

    def __init__(self, p, q=None):
        p = flint.fmpq_poly(p)
        q = flint.fmpq_poly(1) if q is None else flint.fmpq_poly(q)
        if p == 0:
            self.p, self.q = flint.fmpq_poly(0), flint.fmpq_poly(1)
            return
        g = p.gcd(q)
        p, q = p // g, q // g
        lc = q.coeffs()[-1]
        self.p, self.q = p/lc, q/lc

    def __add__(s, o):
        o = o if isinstance(o, RF) else RF(o)
        return RF(s.p*o.q + o.p*s.q, s.q*o.q)
    __radd__ = __add__

    def __neg__(s):
        return RF(-s.p, s.q)

    def __sub__(s, o):
        return s + (-(o if isinstance(o, RF) else RF(o)))

    def __mul__(s, o):
        o = o if isinstance(o, RF) else RF(o)
        return RF(s.p*o.p, s.q*o.q)
    __rmul__ = __mul__

    def __truediv__(s, o):
        o = o if isinstance(o, RF) else RF(o)
        return RF(s.p*o.q, s.q*o.p)

    def iszero(s):
        return s.p == 0

    def at(s, v):
        return Fr(str(s.p(v)))/Fr(str(s.q(v)))

    def __repr__(s):
        return f'({s.p})/({s.q})'


ONE = RF(1)
def sg(a):
    return 1 if a > 0 else -1


def acc(E, w, c):
    if c.iszero():
        return
    if w in E:
        t = E[w] + c
        if t.iszero():
            del E[w]
        else:
            E[w] = t
    else:
        E[w] = c


@lru_cache(maxsize=None)
def qsh(u, v):
    """Z_u * Z_v as {word: int}  -- fact (Q)."""
    if not u:
        return {v: 1}
    if not v:
        return {u: 1}
    out = {}
    def put(w, m):
        out[w] = out.get(w, 0) + m
    for w, m in qsh(u[1:], v).items():
        put((u[0],) + w, m)
    for w, m in qsh(u, v[1:]).items():
        put((v[0],) + w, m)
    ab = sg(u[0])*sg(v[0])*(abs(u[0]) + abs(v[0]))
    for w, m in qsh(u[1:], v[1:]).items():
        put((ab,) + w, m)
    return {w: m for w, m in out.items() if m}


def emul(E1, E2):
    out = {}
    for w1, c1 in E1.items():
        for w2, c2 in E2.items():
            c = c1*c2
            for w, m in qsh(w1, w2).items():
                acc(out, w, c*m)
    return out


def x_at(a, t):
    """x_a(2n + t) as RF (sgn(a)^(2n+t) = sgn(a)^t)."""
    return RF(sg(a)**(t % 2)) / RF((2*X + t)**abs(a))


@lru_cache(maxsize=None)
def Zshift(w, d):
    """Z_w(2n + d) as {word: RF} in Z(2n)  -- fact (P)."""
    if not w or d == 0:
        return {w: ONE}
    out = {}
    if d > 0:
        for ww, c in Zshift(w, d - 1).items():
            acc(out, ww, c)
        xa = x_at(w[0], d)
        for ww, c in Zshift(w[1:], d - 1).items():
            acc(out, ww, xa*c)
    else:
        for ww, c in Zshift(w, d + 1).items():
            acc(out, ww, c)
        xa = x_at(w[0], d + 1)
        for ww, c in Zshift(w[1:], d).items():
            acc(out, ww, -(xa*c))
    return out


@lru_cache(maxsize=None)
def hcoef(w):
    """c_w(N) = sum k sigma^N N^-m Z_v(N-1): tuple of (k Fraction, sigma, m, v)  -- fact (H)."""
    if len(w) == 1:
        return ((Fr(1), 1, 1, ()),) if w[0] == 1 else ((Fr(-1), -1, 1, ()),)
    a = w[0]
    inner = hcoef(w[1:])
    if a == 0:
        return tuple((k, s, m + 1, v) for k, s, m, v in inner)
    return tuple((a*k, a, 1, ((a*s)*m,) + v) for k, s, m, v in inner)


def hcoef_at(w, d):
    """c_w(2n + d) as {word: RF}."""
    out = {}
    for k, s, m, v in hcoef(w):
        pre = RF(k.numerator)/RF(k.denominator)*RF(s**(d % 2))/RF((2*X + d)**m)
        for ww, c in Zshift(v, d - 1).items():
            acc(out, ww, pre*c)
    return out


# ------------------------------------------------------------- exact numeric evaluators (testing only)
@lru_cache(maxsize=None)
def Zv(w, M):
    if not w:
        return Fr(1)
    if M <= 0:
        return Fr(0)
    return Zv(w, M - 1) + Fr(sg(w[0])**M, M**abs(w[0]))*Zv(w[1:], M - 1)


def ev(E, nv):
    return sum((c.at(nv)*Zv(w, 2*nv) for w, c in E.items()), Fr(0))


@lru_cache(maxsize=None)
def hpl_series(w, N):
    """direct Taylor coefficients of H_w from the integral definition (list, index = power)."""
    if len(w) == 1:
        return tuple([Fr(0)] + [Fr(1, m) if w[0] == 1 else Fr((-1)**(m + 1), m) for m in range(1, N + 1)])
    inner = hpl_series(w[1:], N)
    a = w[0]
    out = [Fr(0)]*(N + 1)
    if a == 0:
        for m in range(1, N + 1):
            out[m] = inner[m]/m
    else:
        # f_a(t) H_w(t) = sum_m t^m sum_{i<=m} a^(m-i) inner[i]; integrate
        run = Fr(0)
        for m in range(0, N):
            run = run*a + inner[m]
            out[m + 1] = run/(m + 1)
    return tuple(out)


def Hf(r, m):
    return sum((Fr(1, i**r) for i in range(1, m + 1)), Fr(0))


# ------------------------------------------------------------- data
def load(fname):
    out = []
    for nm, v in json.load(open(G1 + fname)):
        if nm == '1':
            out.append((0, (), Fr(v)))
        elif nm == 'c':
            out.append((1, (), Fr(v)))
        else:
            m = re.fullmatch(r'c\^(-?\d+) H\[([-0-9, ]+)\]', nm)
            assert m, nm
            out.append((int(m.group(1)), tuple(int(t) for t in m.group(2).split(',')), Fr(v)))
    return out


def side_hpl(terms, par):
    """[c^(2n+par)] sum q c^e H_w  for n >= 1 (constant/c terms excluded: they sit at c^0, c^1)."""
    E = {}
    for e, w, qv in terms:
        if not w:
            continue
        for ww, c in hcoef_at(w, par - e).items():
            acc(E, ww, RF(qv.numerator)/RF(qv.denominator)*c)
    return E


def sym_to_RF(expr):
    num, den = sp.fraction(sp.cancel(sp.sympify(expr)))
    nsym = sp.Symbol('n')
    def P(e):
        pl = sp.Poly(sp.expand(e), nsym)
        cs = [Fr(0)]*(pl.degree() + 1)
        for (k,), c in pl.terms():
            c = sp.Rational(c)
            cs[k] = Fr(int(c.p), int(c.q))
        return flint.fmpq_poly([flint.fmpq(c.numerator, c.denominator) for c in cs])
    return RF(P(num), P(den))


def side_mu(formula):
    """formula: {monomial (h1,h2,h3,g1,g2,g3,z2,z3,EG): str coefficient} -> {part: expression}."""
    sub = []
    for r in (1, 2, 3):
        sub.append({(r,): RF(2**(r - 1)), (-r,): RF(2**(r - 1)), (): -(ONE/RF(X**r))})   # H^(r)_{n-1}
    for r in (1, 2, 3):
        sub.append({(r,): ONE, (): -(ONE/RF((2*X)**r))})                                   # H^(r)_{2n-1}
    parts = {'1': {}, 'z2': {}, 'z3': {}}
    for mon, cs in formula.items():
        assert mon[8] == 0
        E = {(): ONE}
        for S, e in zip(sub, mon[:6]):
            for _ in range(e):
                E = emul(E, S)
        key = 'z2' if mon[6] else ('z3' if mon[7] else '1')
        assert mon[6] + mon[7] <= 1
        cf = sym_to_RF(cs)
        for w, c in E.items():
            acc(parts[key], w, cf*c)
    return parts


def mu_value(formula, nv):
    vals = [Hf(1, nv - 1), Hf(2, nv - 1), Hf(3, nv - 1), Hf(1, 2*nv - 1), Hf(2, 2*nv - 1), Hf(3, 2*nv - 1)]
    out = {'1': Fr(0), 'z2': Fr(0), 'z3': Fr(0)}
    nsym = sp.Symbol('n')
    for mon, cs in formula.items():
        cv = sp.Rational(sp.sympify(cs).subs(nsym, nv))
        v = Fr(int(cv.p), int(cv.q))
        for x, e in zip(vals, mon[:6]):
            v *= x**e
        out['z2' if mon[6] else ('z3' if mon[7] else '1')] += v
    return out


def low_orders(terms, part, formula):
    ok = True
    for N in range(-2, 4):
        tot = Fr(0)
        for e, w, qv in terms:
            if not w:
                tot += qv if N == e else 0
            elif N - e >= 1:
                tot += qv*hpl_series(w, N - e)[N - e]
        target = mu_value(formula, N//2)[part] if (N >= 2 and N % 2 == 0) else Fr(0)
        ok = ok and tot == target
    return ok


if __name__ == '__main__':
    random.seed(1)
    # ---- self-tests of (Q), (P), (H), (h) at many n (exact)
    t_ok = True
    words = [(1,), (-1,), (2,), (-3,), (1, -1), (-1, 2), (2, -1, -1), (-1, -1, -1), (3, 1)]
    for nv in (1, 2, 3, 5, 8):
        for u_ in words:
            for v_ in words[:6]:
                if ev({w: RF(m) for w, m in qsh(u_, v_).items()}, nv) != Zv(u_, 2*nv)*Zv(v_, 2*nv):
                    t_ok = False
            for d in range(-2, 4):
                if 2*nv + d >= 0 and ev(Zshift(u_, d), nv) != Zv(u_, 2*nv + d):
                    t_ok = False
    for w in [(1,), (-1,), (0, 1), (1, -1), (-1, 0, 1), (0, -1, 1, -1), (1, 1, 0, -1), (-1, -1, 0, 1)]:
        ser = hpl_series(w, 40)
        for nv in (1, 2, 5, 9, 15):
            for d in range(-1, 4):
                if ev(hcoef_at(w, d), nv) != ser[2*nv + d]:
                    t_ok = False
    print('self-tests of quasi-shuffle, peeling, HPL-coefficient and harmonic-number rules (exact, several n):', t_ok)
    # also check hpl_series against mpmath for a few classical functions
    import mpmath as mp
    mp.mp.dps = 30
    x0 = mp.mpf(1)/3
    def val(w, N=200):
        s = hpl_series(w, N)
        return sum(mp.mpf(c.numerator)/c.denominator*x0**k for k, c in enumerate(s))
    print('  H_{0,0,1}(1/3) - Li3(1/3)       :', mp.nstr(val((0, 0, 1)) - mp.polylog(3, x0), 3))
    print('  H_{0,-1}(1/3) + Li2(-1/3)       :', mp.nstr(val((0, -1)) + mp.polylog(2, -x0), 3))
    print('  H_{1,1}(1/3) - log(1-x)^2/2     :', mp.nstr(val((1, 1)) - mp.log(1 - x0)**2/2, 3))

    mu = pickle.load(open('v_mu_formula.pkl', 'rb'))
    allok = True
    for k, files in ((3, {'1': 'hpl_fit_k3_(0,0,0).json', 'z2': 'hpl_fit_k3_(1,0,0).json', 'z3': 'hpl_fit_k3_(0,1,0).json'}),
                     (2, {'1': 'hpl_fit_k2_(0,0,0).json', 'z2': 'hpl_fit_k2_(1,0,0).json'})):
        MU = side_mu(mu[k])
        if k == 2:
            print(f'B{k}: mu_(2,n) zeta3 part empty: {len(MU["z3"]) == 0}')
        for part, fn in files.items():
            terms = load(fn)
            ev_ = side_hpl(terms, 0)
            od_ = side_hpl(terms, 1)
            dif = dict(ev_)
            for w, c in MU[part].items():
                acc(dif, w, -c)
            lo = low_orders(terms, part, mu[k])
            # spot check the canonical forms against direct numbers
            sp_ok = all(ev(ev_, nv) == sum((qv*hpl_series(w, 2*nv + 3)[2*nv - e] for e, w, qv in terms if w), Fr(0))
                        and ev(MU[part], nv) == mu_value(mu[k], nv)[part] for nv in (1, 2, 4, 7))
            ok = (not dif) and (not od_) and lo and sp_ok
            allok = allok and ok
            print(f'B{k} part {part}: {len(terms)} terms | canonical words even/mu: {len(ev_)}/{len(MU[part])} | '
                  f'even identity: {not dif} | odd identity: {not od_} | c^-2..c^3 exact: {lo} | spot checks: {sp_ok}')
            sys.stdout.flush()
    print('B3 and B2 hold for all n >= 1 (verifier, independent code):', allok)

    # ---- mutation controls: perturb a random coefficient; identity must fail
    for k, fn in ((3, 'hpl_fit_k3_(0,0,0).json'), (3, 'hpl_fit_k3_(1,0,0).json')):
        terms = load(fn)
        MU = side_mu(mu[k])
        part = '1' if '(0,0,0)' in fn else 'z2'
        fails = 0
        for trial in range(5):
            i = random.randrange(len(terms))
            mt = list(terms)
            e, w, qv = mt[i]
            if not w:
                continue
            mt[i] = (e, w, qv + Fr(1, 7))
            ev_ = side_hpl(mt, 0)
            od_ = side_hpl(mt, 1)
            dif = dict(ev_)
            for ww, c in MU[part].items():
                acc(dif, ww, -c)
            fails += int(bool(dif) or bool(od_))
        print(f'mutation control on {fn}: identity fails in {fails} of the perturbed trials (must equal number of trials run)')
    # ---- weight-4 Taylor coefficients of W4(c): canonical form (for the w_n formula)
    R3 = load('hpl_fit_k3_(0,0,0).json'); PSI = load('hpl_fit_k3_(1,0,0).json'); MPHI = load('hpl_fit_k3_(0,1,0).json')
    W4 = {}
    for part, terms, wt in (('1', R3, 4), ('z2', PSI, 2), ('z3', MPHI, 1)):
        W4[part] = side_hpl([t for t in terms if len(t[1]) == wt], 0)
        W4odd = side_hpl([t for t in terms if len(t[1]) == wt], 1)
        print(f'W4 {part}-part: even canonical words {len(W4[part])}, odd part vanishes identically: {not W4odd}')
    pickle.dump({p: {w: (str(c.p), str(c.q)) for w, c in E.items()} for p, E in W4.items()}, open('v_w4_canonical.pkl', 'wb'))
