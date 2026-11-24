# Exact symbolic summation of rational functions over arithmetic progressions,
# reduced to harmonic numbers H^{(s)}_{a t} with t symbolic (n = 4t+1 or 4t+3).
import sympy as sp
t, i = sp.symbols('t i', integer=True, nonnegative=True)
LN2, GAM = sp.Symbol('ln2'), sp.Symbol('gammaE')
ZETA = {s: sp.Symbol('zeta%d' % s) for s in range(2, 8)}
_H = {}
def Hsym(s, a):
    if a == 0: return 0
    key = (s, a)
    if key not in _H: _H[key] = sp.Symbol('H%d_%dt' % (s, a))
    return _H[key]

def H_shift(s, a, k):
    """H^{(s)}_{a t + k} in terms of H^{(s)}_{a t} plus rational terms (k integer)."""
    base = Hsym(s, a)
    if k >= 0: return base + sum(sp.Rational(1)/(a*t+j)**s for j in range(1, k+1))
    return base - sum(sp.Rational(1)/(a*t-j)**s for j in range(0, -k))

def Zred(s, a, g):
    """sum_{i>=0} 1/(i + a t + g)^s for s>=2 ; for s=1 returns -psi(a t + g) (to be combined)."""
    g = sp.Rational(g)
    if g.q == 1:
        N_k = int(g) - 1                    # zeta(s, a t + k) with k = g: = zeta(s) - H_{a t + k - 1}
        if s == 1: return GAM - H_shift(1, a, N_k)
        return ZETA[s] - H_shift(s, a, N_k)
    assert g.q == 2
    k = int(g - sp.Rational(1, 2))          # M = a t + k
    if s == 1: return GAM + 2*LN2 - 2*H_shift(1, 2*a, 2*k) + H_shift(1, a, k)
    return (2**s-1)*ZETA[s] - 2**s*H_shift(s, 2*a, 2*k) + H_shift(s, a, k)

def _terms(f):
    f = sp.apart(sp.together(f), i)
    out = []
    for term in sp.Add.make_args(f):
        num, den = sp.fraction(sp.factor(term))
        # den = const * (p i + q)^s
        c = num; s = 0; lin = None
        for fac, e in sp.factor_list(den)[1]:
            if fac.has(i):
                assert lin is None; lin, s = fac, e
            else:
                c = c / fac**e
        c = c / sp.factor_list(den)[0]
        if lin is None:
            raise ValueError('polynomial part %s' % term)
        p = sp.Poly(lin, i).coeffs()[0]; q = lin - p*i
        a_t = sp.expand(q/p)
        out.append((sp.simplify(c/p**s), s, a_t))
    return out

def _split(a_t):
    a_t = sp.expand(a_t)
    a = sp.Poly(a_t, t).coeff_monomial(t) if a_t.has(t) else 0
    g = a_t - a*t
    return int(a), sp.Rational(g)

def sum_inf(f, i0=0):
    """sum_{i>=i0} f(i)."""
    f = f.subs(i, i + i0)
    res = 0; s1 = 0
    for c, s, a_t in _terms(f):
        a, g = _split(a_t)
        assert a >= 0 and (a > 0 or g > 0), (a, g, f)
        res += c*Zred(s, a, g)
        if s == 1: s1 += c
    assert sp.simplify(s1) == 0, s1
    return sp.expand(res)

def sum_fin(f, i0, i1):
    """sum_{i=i0}^{i1} f(i), i1 may depend on t (i1 = b t + const)."""
    f = f.subs(i, i + i0); U = sp.expand(i1 - i0)
    res = 0
    for c, s, a_t in _terms(f):
        a, g = _split(a_t)
        if a < 0:   # reverse the range: i -> U - j
            c = c*(-1)**s; a_t = sp.expand(-U - a_t); a, g = _split(a_t)
        assert a >= 0 and (a > 0 or g > 0), (a, g)
        a2, g2 = _split(a_t + U + 1)
        res += c*(Zred(s, a, g) - Zred(s, a2, g2))
    return sp.expand(res)
