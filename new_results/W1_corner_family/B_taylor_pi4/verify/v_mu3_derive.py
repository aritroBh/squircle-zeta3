"""v_mu3_derive.py -- verifier's independent symbolic derivation of mu_{k,n}, k <= 3 (G1 Theorem A).

Different route from the solver's mu3_symbolic.py:
  * the u-expansion of the integrand is done by an explicit binomial series
        (1 + y W)^{1-2n} = sum_j binom(1-2n, j) y^j W^j,   W = X - 1 = sum_m u^m l^m/m!,
    times the truncated exponential exp(u A)  (no sympy series of a logarithm);
  * Beta derivatives are obtained by letting sympy differentiate exp(loggamma(p)+loggamma(q)-loggamma(p+q))
    (no hand-written cumulant/moment recursion), then polygamma(m, n+j) and polygamma(m, 2n+j) are rewritten
    with the shift rule and the integer values psi(N) = -gamma + H_{N-1}, psi^(m)(N) = (-1)^(m+1) m!(zeta(m+1) - H^(m+1)_{N-1}).
Integrand (G1 Theorem A, y-form, re-derived by the verifier by hand):
    J_n(u) = 1/2 int_0^1 y^((n-1)(1-u)) (1-y)^(n-1+nu) [1 + y(((1-y)/y)^u - 1)]^(1-2n) dy,
    M(u;alpha) = - sum_n C(2n,n)/(2n-1) c^{2n} J_n(u),  C(2n,n) B(n,n)/2 = 1/n.
Output: mu_{k,n} as a polynomial in h1,h2,h3 (= H^(r)_{n-1}), g1,g2,g3 (= H^(r)_{2n-1}), z2, z3, EulerGamma,
with coefficients in Q(n). Compared (a) symbolically with the formula typed from the solver's result.md, and
(b) exactly with G1's stored table flat_exact_N320_K3.json for every n <= 320.
"""
import json, sys, pickle
from fractions import Fraction as Fr
import sympy as sp

n, u, y, a, b = sp.symbols('n u y a b')
p, q = sp.symbols('p q')
h1, h2, h3, g1, g2, g3, z2, z3, eg = sp.symbols('h1 h2 h3 g1 g2 g3 z2 z3 EG')
KMAX = 3

# ---- integrand expansion
l = b - a
A = -(n - 1)*a + n*b
W = sum(u**m*l**m/sp.factorial(m) for m in range(1, KMAX + 1))
binom = lambda j: sp.Mul(*[(1 - 2*n - t) for t in range(j)])/sp.factorial(j)
base = sum(binom(j)*y**j*W**j for j in range(0, KMAX + 1))
expA = sum((u*A)**m/sp.factorial(m) for m in range(0, KMAX + 1))
E = sp.expand(base*expA)
coef_u = {k: sp.expand(E.coeff(u, k)) for k in range(KMAX + 1)}

# ---- Beta derivative ratios via sympy's own differentiation of loggamma
Bf = sp.exp(sp.loggamma(p) + sp.loggamma(q) - sp.loggamma(p + q))
_dr = {}
def Dratio(i, ll):
    if (i, ll) not in _dr:
        d = Bf
        if i:
            d = sp.diff(d, p, i)
        if ll:
            d = sp.diff(d, q, ll)
        _dr[(i, ll)] = sp.expand(sp.simplify(d/Bf))
    return _dr[(i, ll)]

def psi_int(m, Hsym):
    """psi^(m)(N) for an integer N, with Hsym[r] = H^(r)_{N-1}."""
    if m == 0:
        return -eg + Hsym[1]
    zeta = {1: z2, 2: z3}[m]
    return (-1)**(m + 1)*sp.factorial(m)*(zeta - Hsym[m + 1])

Hn = {1: h1, 2: h2, 3: h3}
H2n = {1: g1, 2: g2, 3: g3}

def psi_shift(m, which, j):
    """psi^(m)(N0 + j), N0 = n or 2n: psi^(m)(x+1) = psi^(m)(x) + (-1)^m m!/x^(m+1)."""
    N0 = n if which == 'n' else 2*n
    Hs = Hn if which == 'n' else H2n
    val = psi_int(m, Hs)
    for t in range(j):
        val += (-1)**m*sp.factorial(m)/(N0 + t)**(m + 1)
    return val

def subs_psi(expr, j):
    reps = {}
    for pg in expr.atoms(sp.polygamma):
        m, arg = pg.args
        m = int(m)
        if arg == p:
            reps[pg] = psi_shift(m, 'n', j)
        elif arg == q:
            reps[pg] = psi_shift(m, 'n', 0)
        elif sp.simplify(arg - (p + q)) == 0:
            reps[pg] = psi_shift(m, '2n', j)
        else:
            raise ValueError(arg)
    return expr.xreplace(reps)

def mu(k):
    P = sp.Poly(coef_u[k], y, a, b)
    tot = 0
    for (j, i, ll), co in P.terms():
        poch = sp.Mul(*[(n + t)/(2*n + t) for t in range(j)])   # B(n+j,n)/B(n,n)
        tot += co*poch*subs_psi(Dratio(i, ll), j)
    tot = -tot/(n*(2*n - 1))
    return sp.Poly(sp.expand(tot), h1, h2, h3, g1, g2, g3, z2, z3, eg)

def clean(Pl):
    """dict monomial -> cancelled rational function of n"""
    out = {}
    for mon, cf in Pl.terms():
        cf = sp.factor(sp.cancel(cf))
        if cf != 0:
            out[mon] = cf
    return out

# ---- solver's formula as typed in result.md section 2 (verifier's transcription)
D = g1 - h1
solver_mu3 = (-z3/(4*n**2 - 1)
              + z2*((5*n**2 + 5*n + 1)/(2*n*(n + 1)*(2*n + 1)**2) - D/(4*n**2 - 1))
              + D**3/(6*n*(2*n - 1)) - D**2/(4*n**2) + D*g2/(2*n*(2*n - 1)) - D*h2/(2*n*(4*n**2 - 1))
              + (4*n**2 - 2*n - 1)*D/(4*n**3*(2*n + 1)) + g3/(3*n*(2*n - 1)) + (n - 1)*h3/(3*n*(4*n**2 - 1))
              - g2/(4*n**2) - (3*n + 1)*(2*n**2 - 1)*h2/(4*n**2*(n + 1)*(2*n + 1)**2)
              - (16*n**5 - 12*n**3 + 4*n**2 + 5*n + 1)/(8*n**4*(n + 1)*(2*n + 1)**2))
# G1's mu_{2,n} (G1 result.md section 2)
G1_mu2 = (z2/(4*n**2 - 1) - D**2/(2*n*(2*n - 1)) + D/(2*n**2) + h2/(2*n*(4*n**2 - 1)) - g2/(2*n*(2*n - 1))
          - (4*n**2 - 2*n - 1)/(4*n**3*(2*n + 1)))

def Hf(r, m):
    return sum((Fr(1, i**r) for i in range(1, m + 1)), Fr(0))

if __name__ == '__main__':
    NCHK = int(sys.argv[1]) if len(sys.argv) > 1 else 320
    res = {}
    for k in range(KMAX + 1):
        res[k] = clean(mu(k))
        has_gamma = any(mon[8] for mon in res[k])
        print(f'mu_{{{k},n}}: {len(res[k])} monomials; EulerGamma present: {has_gamma}')
    # symbolic comparison with the typed formulas
    gens = (h1, h2, h3, g1, g2, g3, z2, z3, eg)
    for k, ref in ((3, solver_mu3), (2, G1_mu2)):
        mine = sum(cf*sp.Mul(*[v**e for v, e in zip(gens, mon)]) for mon, cf in res[k].items())
        dif = sp.Poly(sp.expand(sp.together(mine - ref)*1), *gens)
        dd = [sp.cancel(cf) for _, cf in dif.terms()]
        print(f'k={k}: verifier derivation minus typed formula == 0 identically in n: {all(x == 0 for x in dd)}')
    print('mu_{0,n} =', res[0])
    print('mu_{1,n} =', res[1])
    # exact comparison with G1 table
    tab = json.load(open('research/10_frontier/G1_zeta3_any_angle/flat_exact_N320_K3.json'))
    keymap = {'(0, 0, 0)': (0, 0), '(1, 0, 0)': (1, 0), '(0, 1, 0)': (0, 1)}
    lamb = {}
    for k in range(KMAX + 1):
        lamb[k] = [(mon, (lambda cf_: (lambda v: cf_.subs(n, v)))(cf)) for mon, cf in res[k].items()]
    for k in range(KMAX + 1):
        ok = True
        for nv in range(1, NCHK + 1):
            vals = [Hf(1, nv - 1), Hf(2, nv - 1), Hf(3, nv - 1), Hf(1, 2*nv - 1), Hf(2, 2*nv - 1), Hf(3, 2*nv - 1)]
            out = {(0, 0): Fr(0), (1, 0): Fr(0), (0, 1): Fr(0)}
            for mon, f in lamb[k]:
                cv = sp.Rational(f(sp.Integer(nv)))
                v = Fr(int(cv.p), int(cv.q))
                for x, e in zip(vals, mon[:6]):
                    v *= x**e
                out[(mon[6], mon[7])] += v
            g = tab[str(k)][str(nv)]
            gv = {keymap[kk]: Fr(vv) for kk, vv in g.items()}
            for kk in out:
                if out[kk] != gv.get(kk, Fr(0)):
                    ok = False
            if set(gv) - set(out):
                ok = False
        print(f'k={k}: verifier formula == G1 exact table for all n <= {NCHK}: {ok}')
        sys.stdout.flush()
    pickle.dump({k: {mon: str(cf) for mon, cf in res[k].items()} for k in res}, open('v_mu_formula.pkl', 'wb'))
