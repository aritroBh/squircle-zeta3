"""Referee's own exact derivation of the d-expansion of the arc length of one edge.
Radial projection x_p(t) = y(t) F^{-1/p}, F = sum_j a_j(t)^p, a_j affine in t, sum_j a_j = 1 on the edge.
X_j stands for log a_j(t). Everything is a truncated series in d (order 3)."""
import sympy as sp, pickle, sys
t, d, L = sp.symbols('t d L', positive=True)
al, be, B1, C1 = sp.symbols('alpha beta B1 C1')
x0, y0, x1, y1 = sp.symbols('x0 y0 x1 y1', real=True)
X0, X1, Xf1, Xf2 = sp.symbols('X0 X1 Xf1 Xf2')
N = 4
def trunc(e):
    e = sp.expand(e)
    return [e.coeff(d, k) for k in range(N)]
def mul(a, b):
    return [sp.expand(sum(a[k]*b[n-k] for k in range(n+1))) for n in range(N)]
def series_exp0(f):  # f[0] = 0
    g = [sp.Integer(1)] + [0]*(N-1)
    for n in range(1, N):
        g[n] = sp.expand(sum(k*f[k]*g[n-k] for k in range(1, n+1))/n)
    return g
def series_log1(f):  # f[0] = 1
    g = [sp.Integer(0)] + [0]*(N-1)
    for n in range(1, N):
        g[n] = sp.expand(f[n] - sum(k*g[k]*f[n-k] for k in range(1, n))/n)
    return g
def series_inv(b):  # b[0] = c
    q = [0]*N
    for n in range(N):
        q[n] = sp.expand((( 1 if n == 0 else 0) - sum(q[k]*b[n-k] for k in range(n)))/b[0])
    return q
def series_sqrt1(f):  # f[0]=1
    g = [sp.Integer(1)] + [0]*(N-1)
    for n in range(1, N):
        g[n] = sp.expand((f[n] - sum(g[k]*g[n-k] for k in range(1, n)))/2)
    return g

def build(terms):
    """terms: list of (a(t), X symbol). returns [E0..E3] (coefficients of d^k of |x_p'|)."""
    F = [0]*N; G = [0]*N
    for a, X in terms:
        e = [X**k/sp.factorial(k) for k in range(N)]
        ap = sp.diff(a, t)
        for k in range(N):
            F[k] += a*e[k]; G[k] += ap*e[k]
    F = [sp.expand(f) for f in F]; G = [sp.expand(g) for g in G]
    F[0] = sp.simplify(F[0])
    assert sp.simplify(F[0]-1) == 0, F[0]
    F[0] = sp.Integer(1)
    Phi = mul(G, series_inv(F))
    yv = (x0+t*(x1-x0), y0+t*(y1-y0)); yp = (x1-x0, y1-y0)
    yyp = sp.expand(yv[0]*yp[0]+yv[1]*yp[1]); yy = sp.expand(yv[0]**2+yv[1]**2)
    Phi2 = mul(Phi, Phi)
    Q = [sp.expand(-2*Phi[k]*yyp + Phi2[k]*yy) for k in range(N)]  # plus |y'|^2 = L^2
    q = [sp.Integer(1)] + [sp.expand(Q[k]/L**2) for k in range(1, N)]
    assert Q[0] == 0 or sp.simplify(Q[0]) == 0  # Phi0 = sum a_j' = 0
    sq = series_sqrt1(q)
    lf = series_log1(F)
    inv1pd = [(-1)**k for k in range(N)]
    E = series_exp0([-c for c in mul(lf, inv1pd)])
    return [sp.expand(L*c) for c in mul(E, sq)]

if __name__ == '__main__':
    terms = [(al*t, X0), (be*(1-t), X1), (B1+C1*t, Xf1), (1-al*t-be*(1-t)-B1-C1*t, Xf2)]
    E = build(terms)
    pickle.dump(E, open('edge_E_general.pkl', 'wb'))
    E3 = sp.Poly(E[3], X0, X1, Xf1, Xf2)
    g21 = E3.coeff_monomial(X0**2*X1); g12 = E3.coeff_monomial(X0*X1**2)
    print('free symbols of g21:', sorted(map(str, g21.free_symbols)))
    print('free symbols of g12:', sorted(map(str, g12.free_symbols)))
    Lsq = (x1-x0)**2+(y1-y0)**2
    I = sp.integrate(g21+g12, (t, 0, 1))
    Zder = sp.expand(12*I)
    tau_dot_V0 = ((x1-x0)*x0+(y1-y0)*y0)/L
    L0 = -tau_dot_V0; L1 = L-L0; h2 = x0**2+y0**2-L0**2
    Zpap = 6*al*be*h2/L*(3*(al*L0+be*L1)/L-2)
    diff = sp.expand((Zder-Zpap)*L**5)
    diff = sp.expand(diff.subs(L**2, Lsq))  # reduce even powers
    for k in range(6, 1, -1):
        diff = sp.expand(diff.subs(L**k, Lsq**(k//2)*L**(k % 2)))
    print('Z_derived - Z_paper (times L^5, with L^2=|V1-V0|^2):', sp.simplify(diff))
    print('Z_derived =', sp.factor(Zder))
    # first-order coefficient: check it is a total derivative
    E1 = sp.Poly(E[1], X0, X1, Xf1, Xf2)
    print('E1 monomials:', E1.monoms())
