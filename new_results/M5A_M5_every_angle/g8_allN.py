"""g8_allN.py -- the exact Z-sum identity behind Theorems 20 and 21 (sections g8:sec:M2alpha, g8:sec:M4alpha).

What it does (usage: python g8_allN.py K, with K in {2,3,4}).
  (1) mu side.  From Definition g6:def:formal (sec_g6_flat.tex) it computes the formal coefficient P_{K,n} for
      SYMBOLIC n: the polynomials e^{(K)}_{j,i,m}(n) of Lemma g6:lem:En (series in u, ring Q[n,y,a,b]), the
      polynomials D^#_{i,m} in formal harmonic sums h^(r)_l and formal zeta values, and the substitution of the
      algebra D of Lemma g6:lem:A (harmonic numbers -> combinations of alternating Z-sums Z_w(2n) with
      coefficients in Q(n)).  Result, for each monomial of formal zeta values: a dict {word: element of Q(n)}.
  (2) HPL side.  For every term q c^e H_w (and the 'He' terms q c^e (H_w(c) + s H_{-w}(c))) it computes
      [c^N] with N = 2n and N = 2n+1 as an element of the same form (Lemma g6:lem:dict, window shifts).
  (3) It compares the two dicts WORD BY WORD (equality in Q(n), nothing else), and checks the finitely many
      remaining coefficients N <= 1 in exact rational arithmetic.
  (4) As a test of (1) only, it evaluates the Z-sum form of P_{K,n} at n = 1..8 exactly and compares with
      the independent routine g6_common.py (numeric n).
Everything is exact (sympy field Q(n), python Fractions).  Written for this supplement; the Z-sum library is
g8_zsum.py (copy of the notes' zsum.py).
"""
import sys, json, re, itertools, time
from fractions import Fraction as Fr
from math import factorial
from sympy import QQ, ring
sys.setrecursionlimit(20000)
import g8_zsum as zs
from g8_zsum import K, nn, ONE, add_into, lin, mul, Z_val

# ------------------------------------------------------------------ (1) the mu side

def E_coeff_polys(k):
    """[u^k] E_n as dict {(j,i,m): polynomial in n (sympy PolyElement in ring Q[n])}."""
    R, u, y, a, b, n = ring("u,y,a,b,n", QQ)
    def trunc(p):
        return R.from_dict({m: c for m, c in p.items() if m[0] <= k})
    ell = b - a
    A = -(n - 1)*a + n*b
    X = R.zero; term = R.one
    for r in range(1, k + 1):
        term = trunc(term*u*ell)*QQ(1, r)
        X = X + y*term
    L = R.zero; Xp = R.one
    for m in range(1, k + 1):
        Xp = trunc(Xp*X)
        L = L + Xp*QQ((-1)**(m - 1), m)
    logE = u*A + (1 - 2*n)*L
    E = R.one; P = R.one
    for t in range(1, k + 1):
        P = trunc(P*logE)*QQ(1, t)
        E = E + P
    out = {}
    for (du, dy, da, db, dn), c in E.items():
        if du == k:
            out.setdefault((dy, da, db), {})[dn] = Fr(int(c.numerator), int(c.denominator))
    return {key: sum((K(cf.numerator)/K(cf.denominator)*nn**d for d, cf in poly.items()), K(0))
            for key, poly in out.items()}


def exp_Delta(k):
    """exp(Delta_*) truncated at total degree k: dict {(i, m): {(hmon, zmon): Fraction}}.
    hmon = exponent tuple over h1_1..h1_k,h2_1..h2_k,h3_1..h3_k ; zmon = exponent tuple over zeta_2..zeta_k."""
    names = ["x", "xp"] + [f"h{l}_{r}" for l in (1, 2, 3) for r in range(1, k + 1)] + [f"z{w}" for w in range(2, k + 1)]
    R = ring(",".join(names), QQ)
    Rr, gens = R[0], R[1:]
    x, xp = gens[0], gens[1]
    H = {(l, r): gens[2 + (l - 1)*k + (r - 1)] for l in (1, 2, 3) for r in range(1, k + 1)}
    Zs = {w: gens[2 + 3*k + (w - 2)] for w in range(2, k + 1)}
    def trunc(p):
        return Rr.from_dict({m: c for m, c in p.items() if m[0] + m[1] <= k})
    s = x + xp
    Delta = Rr.zero
    for r in range(1, k + 1):
        Delta = Delta + QQ((-1)**(r - 1), r)*(H[(1, r)]*x**r + H[(2, r)]*xp**r - H[(3, r)]*s**r)
    for w in range(2, k + 1):
        Delta = Delta + QQ((-1)**(w + 1), w)*Zs[w]*(s**w - x**w - xp**w)
    E = Rr.one; P = Rr.one
    for t in range(1, k + 1):
        P = trunc(P*Delta)*QQ(1, t)
        E = E + P
    out = {}
    for mon, c in E.items():
        i, m = mon[0], mon[1]
        hm = mon[2:2 + 3*k]; zm = mon[2 + 3*k:]
        cf = Fr(int(c.numerator), int(c.denominator))*factorial(i)*factorial(m)
        out.setdefault((i, m), {})[(tuple(hm), tuple(zm))] = cf
    return out


def chi(l, r, j):
    """element of D (as Z-dict at N = 2n) for H^(r)_{n+j-1} (l=1), H^(r)_{n-1} (l=2), H^(r)_{2n+j-1} (l=3)."""
    if l == 2:
        base = zs.harm_n_minus_1(r)
        return dict(base)
    if l == 1:
        E = dict(zs.harm_n_minus_1(r))
        for t in range(j):
            add_into(E, (), K(1)/(nn + t)**r)
        return E
    E = dict(zs.harm_2n_minus_1(r))
    for t in range(j):
        add_into(E, (), K(1)/(2*nn + t)**r)
    return E


def rho(j):
    out = K(1)
    for t in range(j):
        out = out*(nn + t)/(2*nn + t)
    return out


def mu_side(k):
    """{zmon: {word: K}} = formal coefficients P_{k,n} as Z-sum combinations at 2n."""
    Epoly = E_coeff_polys(k)
    DE = exp_Delta(k)
    pref = -K(1)/(nn*(2*nn - 1))
    res = {}
    cache = {}
    for (j, i, m), epol in Epoly.items():
        assert i + m == k
        for (hm, zm), cf in DE[(i, m)].items():
            key = (j, hm)
            if key not in cache:
                E = {(): K(1)}
                for idx, ex in enumerate(hm):
                    if ex:
                        l = idx//k + 1; r = idx % k + 1
                        c_el = chi(l, r, j)
                        for _ in range(ex):
                            E = mul(E, c_el)
                cache[key] = E
            E = cache[key]
            f = pref*epol*rho(j)*K(cf.numerator)/K(cf.denominator)
            tgt = res.setdefault(zm, {})
            for w, c in E.items():
                add_into(tgt, w, f*c)
    return res


def zmon_name(zm):
    parts = [f"z{w}" + (f"^{e}" if e > 1 else "") for w, e in zip(range(2, 2 + len(zm)), zm) if e]
    return "*".join(parts) if parts else "1"

# ------------------------------------------------------------------ (2) the HPL side

def nz(w):
    return sum(1 for a in w if a != 0)

def neg(w):
    return tuple(-a for a in w)

def parse_terms(fname):
    """terms (kind, w, e, q): kind 'plain' = q c^e H_w(c); 'even' = q c^e (H_w(c) + (-1)^(e+nz) H_{-w}(c))."""
    out = []
    for nm, v in json.load(open(fname)):
        if nm in ('1', 'c'):
            out.append(('const', (), 0 if nm == '1' else 1, Fr(v))); continue
        m = re.match(r'c\^(-?\d+) (He?)\[(.*)\]', nm)
        if m is None:
            m = re.match(r'c\^(-?\d+)$', nm)
            raise ValueError(nm)
        e = int(m.group(1)); w = tuple(int(t) for t in m.group(3).split(','))
        out.append(('even' if m.group(2) == 'He' else 'plain', w, e, Fr(v)))
    return out


def hpl_side(terms, parity):
    E = {}
    for kind, w, e, q in terms:
        if kind == 'const':
            continue
        assert e <= 1
        qK = K(q.numerator)/K(q.denominator)
        for ww, c in zs.hpl_coeff_at(w, parity - e).items():
            add_into(E, ww, qK*c)
        if kind == 'even':
            s = 1 if (e + nz(w)) % 2 == 0 else -1
            for ww, c in zs.hpl_coeff_at(neg(w), parity - e).items():
                add_into(E, ww, qK*s*c)
    return E


def hpl_coeffs(w, order):
    """exact [x^N] H_w(x), N = 0..order (own recursion from the definition H_{a,w'} = int f_a H_w')."""
    if len(w) == 1:
        a = w[0]
        return [Fr(0)] + [Fr(1, m) if a == 1 else Fr((-1)**(m + 1), m) for m in range(1, order + 1)]
    inner = hpl_coeffs(w[1:], order)
    a = w[0]
    out = [Fr(0)]*(order + 1)
    if a == 0:
        for m in range(1, order + 1):
            out[m] = inner[m]/m
    else:
        acc = Fr(0)   # (1/(1 -/+ x)) * inner : partial sums with sign a
        part = []
        for m in range(order + 1):
            acc = acc*a + inner[m]
            part.append(acc)       # coefficient of x^m in inner / (1 - a x)
        for m in range(1, order + 1):
            out[m] = part[m - 1]/m
    return out


def low_order(terms, target, lo=-4, hi=1):
    """[c^N] of the term list for lo <= N <= hi, exact; compare with target(N) (a Fraction)."""
    bad = []
    for N in range(lo, hi + 1):
        tot = Fr(0)
        for kind, w, e, q in terms:
            if kind == 'const':
                if N == e: tot += q
                continue
            idx = N - e
            if idx >= 1:
                c1 = hpl_coeffs(w, idx)[idx]
                if kind == 'plain':
                    tot += q*c1
                else:
                    s = 1 if (e + nz(w)) % 2 == 0 else -1
                    tot += q*(c1 + s*hpl_coeffs(neg(w), idx)[idx])
        if tot != target(N):
            bad.append((N, tot, target(N)))
    return bad

# ------------------------------------------------------------------ evaluation of Z-sum dicts at integer n

def ev(E, nv):
    tot = Fr(0)
    for w, cf in E.items():
        cv = Fr(int(cf.numer(nv)), 1)/Fr(int(cf.denom(nv)), 1)
        tot += cv*Z_val(w, 2*nv)
    return tot


def run_k(k, parts, artefacts):
    t0 = time.time()
    P = mu_side(k)
    print(f"k={k}: mu side built in {time.time()-t0:.1f}s; formal zeta monomials: " +
          ", ".join(f"{zmon_name(zm)}({len(E)} words)" for zm, E in sorted(P.items())))
    # test of the mu side against the independent routine g6_common (numeric n)
    import g6_common as gc
    okt = True
    for nv in range(1, 9):
        tab = gc.mu_formal(nv, k)[k]
        for zm, E in P.items():
            val = ev(E, nv)
            if val != tab.get(zm, Fr(0)):
                okt = False
                print('   MISMATCH vs g6_common at n =', nv, zmon_name(zm), val, tab.get(zm, Fr(0)))
    print("   mu side equals g6_common.py (independent, numeric n) for n=1..8:", okt)
    allok = True
    for pname, (zm, termfile, extra) in parts.items():
        muE = P.get(zm, {})
        terms = termfile if isinstance(termfile, list) else parse_terms(termfile)
        ev_even = hpl_side(terms, 0)
        for ww, c in extra.items():
            add_into(ev_even, ww, c)
        diff = lin((1, ev_even), (-1, muE))
        odd = hpl_side(terms, 1)
        # target of the low-order check: mu_{k,N/2} part for N>=2 even is covered by the all-n identity; here N <= 1
        lo = low_order(terms, lambda N: Fr(0), lo=-4, hi=1 if not extra else 1)
        # the extra rational terms (Li_2 piece) only live in even N >= 2
        from collections import Counter
        cw = Counter(len(w) for kind, w, e, q in terms if kind != 'const')
        es = sorted({e for kind, w, e, q in terms if kind != 'const'})
        nconst = sum(1 for kind, w, e, q in terms if kind == 'const')
        print(f"   part {pname}: {len(terms)} terms ({nconst} constant); numbers of terms by word length {dict(sorted(cw.items()))}; exponents e in {es}")
        print(f"      HPL side {len(ev_even)} words, mu side {len(muE)} words; depth of the words <= {max((len(w) for w in muE), default=0)}, indices {sorted({a for w in muE for a in w})}")
        print(f"      even identity [c^(2n)] HPL side - mu_(k,n) == 0 as an element of Q(n)-combinations of Z-words, all n>=1: {len(diff) == 0}")
        print(f"      odd identity  [c^(2n+1)] HPL side == 0, all n>=1: {len(odd) == 0}")
        print(f"      [c^N], N=-4..1, exact rational arithmetic, all zero: {len(lo) == 0}")
        if diff: print('      residual:', list(diff.items())[:3])
        if odd: print('      odd residual:', list(odd.items())[:3])
        if lo: print('      low mismatch:', lo)
        allok = allok and not diff and not odd and not lo
        artefacts[pname] = (len(terms), len(ev_even), len(muE))
    return P, allok

# ------------------------------------------------------------------ the three runs

M2_TERMS = None

def m2_terms_from_statement():
    """Y(x) of Theorem 20 of main.tex, expanded in c^e H_w; Y(c)+Y(-c) = sum of 'even' terms with the same coefficients.
    Y = (1-x)[H001 + H01-1 - H1-11] + (1-x^2)/(2x) H101 - (1-x)^2/(2x) H10-1."""
    T = []
    def add(e, w, q): T.append(('even', tuple(w), e, Fr(q)))
    for w, s in (((0, 0, 1), 1), ((0, 1, -1), 1), ((1, -1, 1), -1)):
        add(0, w, s); add(1, w, -s)
    # (1-x^2)/(2x) = 1/(2x) - x/2
    add(-1, (1, 0, 1), Fr(1, 2)); add(1, (1, 0, 1), Fr(-1, 2))
    # -(1-x)^2/(2x) = -1/(2x) + 1 - x/2
    add(-1, (1, 0, -1), Fr(-1, 2)); add(0, (1, 0, -1), 1); add(1, (1, 0, -1), Fr(-1, 2))
    return T


def selftest():
    """the rewriting rules of the HPL side: h_w(2n+delta) (Lemma g6:lem:dict + window shifts) against the exact Taylor
    coefficients computed from the definition of H_w, for every word that occurs in the term lists, delta = -1..4, n = 1..7;
    also the stuffle product and the harmonic-number rewriting."""
    words = set()
    for f in ('g8_terms_k3_part1.json', 'g8_terms_k3_partz2.json', 'g8_terms_k4_part1.json', 'g8_terms_k4_partz2.json'):
        for kind, w, e, q in parse_terms(f):
            if kind != 'const':
                words.add(w); words.add(neg(w))
    for kind, w, e, q in m2_terms_from_statement():
        words.add(w); words.add(neg(w))
    ok = True; cnt = 0
    for w in sorted(words):
        co = hpl_coeffs(w, 2*7 + 4 + 1)
        for nv in range(1, 8):
            for delta in range(-1, 5):
                cnt += 1
                if ev(zs.hpl_coeff_at(w, delta), nv) != co[2*nv + delta]:
                    ok = False; print('FAIL', w, nv, delta)
    for nv in range(1, 8):
        for r in (1, 2, 3, 4):
            Hn = sum((Fr(1, i**r) for i in range(1, nv)), Fr(0)); H2 = sum((Fr(1, i**r) for i in range(1, 2*nv)), Fr(0))
            if ev(zs.harm_n_minus_1(r), nv) != Hn or ev(zs.harm_2n_minus_1(r), nv) != H2:
                ok = False; print('harm FAIL', r, nv)
        for u in [(1,), (-1, 2), (2, -1, 1)]:
            for v in [(2,), (-1, -1), (3, 1)]:
                if ev(mul({u: ONE}, {v: ONE}), nv) != Z_val(u, 2*nv)*Z_val(v, 2*nv):
                    ok = False; print('stuffle FAIL', u, v, nv)
    print(f"selftest: {len(words)} words, {cnt} (word, n, delta) checks of h_w(2n+delta); harmonic numbers; stuffle: all equal = {ok}")


if __name__ == '__main__':
    if sys.argv[1] == 'selftest':
        selftest(); sys.exit(0)
    k = int(sys.argv[1])
    art = {}
    if k == 2:
        # Y(c)+Y(-c) is the 'even' sum.  Add the term -1/2 Li_2(c^2): coefficient -1/(2 n^2) at c^{2n}.
        zm_1 = (0,); zm_z2 = (1,)
        parts = {'1': (zm_1, m2_terms_from_statement(), {(): -K(1)/(2*nn**2)})}
        P, ok = run_k(2, parts, art)
        # zeta(2)-part: [z2] P_{2,n} == 1/(4n^2-1)
        okz = (P[zm_z2] == {(): K(1)/(4*nn**2 - 1)})
        print("   [zeta(2)] P_{2,n} == 1/(4n^2-1) as an element of Q(n):", okz)
        print("THEOREM 20 identity (all n):", ok and okz)
    elif k == 3:
        parts = {'1': ((0, 0), 'g8_terms_k3_part1.json', {}),
                 'z2': ((1, 0), 'g8_terms_k3_partz2.json', {})}
        P, ok = run_k(3, parts, art)
        okz = (P[(0, 1)] == {(): -K(1)/(4*nn**2 - 1)})
        print("   [zeta(3)] P_{3,n} == -1/(4n^2-1):", okz)
        print("R_3 identity (all n):", ok and okz)
    elif k == 4:
        parts = {'1': ((0, 0, 0), 'g8_terms_k4_part1.json', {}),
                 'z2': ((1, 0, 0), 'g8_terms_k4_partz2.json', {})}
        P, ok = run_k(4, parts, art)
        P3 = mu_side(3)
        okl = (P[(0, 1, 0)] == {w: -c for w, c in P3[(1, 0)].items()})
        print("   ladder rung: [zeta(3)] P_{4,n} == -[zeta(2)] P_{3,n} (as elements of Q(n)-Z-word combinations):", okl)
        import sympy as sp
        n_ = sp.Symbol('n')
        for name, zm in (('zeta(4)', (0, 0, 1)), ('zeta(2)^2', (2, 0, 0))):
            E = P[zm]
            print(f"   [{name}] P_(4,n) =", {w: sp.factor(c.as_expr()) for w, c in E.items()})
        print("THEOREM 21 identity (all n, parts 1 and zeta(2)):", ok, "; ladder:", okl)
