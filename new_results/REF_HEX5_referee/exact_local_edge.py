"""Referee's own exact P^{(5)}(1) contribution of a fully local edge (kappa = 0: the far term is constant
on the edge), e.g. every edge of a parallelogram or of an affine-regular hexagon.

Edge y(t) = V0 + t D.  Terms at p = 1+d:  A0 t e^{d(X+m0)},  A1 (1-t) e^{d(Y+m1)},  B e^{d mB},
X = log t, Y = log(1-t), m_j = log(term value / weight).  Boundary point x(t) = y F^{-1/p}, so
|x'| = F^{-1/p} sqrt(L2 - (2/p)(y.D) F'/F + |y|^2 F'^2/(p^2 F^2)),  F' = (1+d)(A0 e^{d(X+m0)} - A1 e^{d(Y+m1)}).
Everything is expanded exactly as a truncated power series in d; the d^k coefficient is a polynomial in
t, X, Y, m0, m1, mB and is integrated exactly with derivatives of the Beta function (sympy).
Returns E_k = k! * int_0^1 [d^k] |x'| dt, k = 0..5.
"""
import sympy as sp
from functools import lru_cache

t, X, Y = sp.symbols('t X Y')
m0, m1, mB = sp.symbols('m0 m1 mB')
NO = 6


def smul(a, b):
    return [sp.expand(sum(a[i] * b[k - i] for i in range(k + 1))) for k in range(NO)]


def sexp_nc(a):
    """exp of a series with a[0] == 0."""
    assert a[0] == 0
    res = [sp.Integer(1)] + [sp.Integer(0)] * (NO - 1)
    term = [sp.Integer(1)] + [sp.Integer(0)] * (NO - 1)
    for k in range(1, NO):
        term = [x / k for x in smul(term, a)]
        res = [sp.expand(res[i] + term[i]) for i in range(NO)]
    return res


def slog1p(u):
    """log(1+u), u[0] == 0."""
    assert u[0] == 0
    res = [sp.Integer(0)] * NO
    pw = [sp.Integer(1)] + [sp.Integer(0)] * (NO - 1)
    for k in range(1, NO):
        pw = smul(pw, u)
        res = [sp.expand(res[i] + sp.Rational((-1) ** (k + 1), k) * pw[i]) for i in range(NO)]
    return res


def ssqrt1p(u):
    assert u[0] == 0
    res = [sp.Integer(1)] + [sp.Integer(0)] * (NO - 1)
    pw = [sp.Integer(1)] + [sp.Integer(0)] * (NO - 1)
    for k in range(1, NO):
        pw = smul(pw, u)
        c = sp.binomial(sp.Rational(1, 2), k)
        res = [sp.expand(res[i] + c * pw[i]) for i in range(NO)]
    return res


def eser(Z):
    return [sp.expand(Z ** k / sp.factorial(k)) for k in range(NO)]


INV1PD = [sp.Integer((-1) ** k) for k in range(NO)]
ONEPD = [sp.Integer(1), sp.Integer(1)] + [sp.Integer(0)] * (NO - 2)


@lru_cache(None)
def beta_moment(n, a, b):
    """int_0^1 t^n log^a t log^b (1-t) dt, exact."""
    x, y = sp.symbols('x y', positive=True)
    Bf = sp.gamma(x) * sp.gamma(y) / sp.gamma(x + y)
    d = Bf
    if a:
        d = sp.diff(d, x, a)
    if b:
        d = sp.diff(d, y, b)
    v = sp.expand_func(d.subs({x: n + 1, y: 1}))
    v = sp.expand(v)
    assert not v.has(sp.EulerGamma), (n, a, b)
    return v


def edge_exact(A0, A1, B, L2, q0, V02):
    e0 = eser(X + m0)
    e1 = eser(Y + m1)
    eB = eser(mB)
    F = [sp.expand(A0 * t * e0[k] + A1 * (1 - t) * e1[k] + B * eB[k]) for k in range(NO)]
    assert sp.simplify(F[0] - 1) == 0, F[0]
    F[0] = sp.Integer(1)
    Fp = smul(ONEPD, [sp.expand(A0 * e0[k] - A1 * e1[k]) for k in range(NO)])
    lF = slog1p([sp.Integer(0)] + F[1:])
    rho = sexp_nc([-c for c in smul(lF, INV1PD)])          # F^{-1/p}
    invF = sexp_nc([-c for c in lF])                          # 1/F
    g = smul(smul(Fp, invF), INV1PD)                          # F'/(pF)
    yD = q0 + t * L2
    yy = V02 + 2 * t * q0 + t ** 2 * L2
    inner = [sp.expand((-2 * yD * g[k] + yy * smul(g, g)[k]) / L2) for k in range(NO)]
    assert inner[0] == 0
    sq = ssqrt1p(inner)
    integrand = smul(rho, sq)          # times L = sqrt(L2)
    out = []
    for k in range(NO):
        P = sp.Poly(integrand[k], t, X, Y)
        tot = sp.Integer(0)
        for (n, a, b), c in P.terms():
            tot += c * beta_moment(n, a, b)
        out.append(sp.expand(sp.factorial(k) * sp.sqrt(L2) * tot))
    return out


if __name__ == '__main__':
    import time, pickle
    t0 = time.time()
    # regular hexagon edge: corners (0,1/sqrt3) -> (-1/2, 1/(2 sqrt3)); A0 = A1 = B = 1/2, weights 1
    s3 = sp.sqrt(3)
    V0 = sp.Matrix([0, 1 / s3]); V1 = sp.Matrix([-sp.Rational(1, 2), 1 / (2 * s3)])
    D = V1 - V0
    E = edge_exact(sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 2), sp.simplify(D.dot(D)),
                   sp.simplify(V0.dot(D)), sp.simplify(V0.dot(V0)))
    lam = sp.log(2)
    sub = {m0: -lam, m1: -lam, mB: -lam}
    P = [sp.expand(6 * e.subs(sub)) for e in E]
    print('time %.1f s' % (time.time() - t0))
    for k in range(NO):
        print('P^(%d)(1) =' % k, sp.N(P[k], 40))
    pi, z3, z5 = sp.pi, sp.zeta(3), sp.zeta(5)
    printed = (s3 / 128) * (30720 * lam - 4920 * pi ** 2 - 61440 * lam ** 2 + 6240 * pi ** 2 * lam + 30720 * lam ** 3
                            + 448 * pi ** 4 + 5760 * z3 * lam - 1920 * pi ** 2 * lam ** 2 - 5120 * lam ** 4 + 11340 * z5
                            + 5580 * (pi ** 2 / 6) * z3 - 123 * pi ** 4 * lam - 1440 * z3 * lam ** 2 + 160 * pi ** 2 * lam ** 3
                            + 256 * lam ** 5)
    diff = sp.simplify(sp.expand(P[5] - printed))
    print('regular hexagon: exact P5 - printed closed form =', diff)
    printed3 = s3 / 8 * (16 * lam ** 3 + 3 * pi ** 2 * lam - 96 * lam ** 2 + 96 * lam - 9 * z3 - 6 * pi ** 2)
    print('regular hexagon: exact P3 - printed closed form =', sp.simplify(sp.expand(P[3] - printed3)))
    print('P1 - 2 sqrt3 log2 =', sp.simplify(P[1] - 2 * s3 * lam))
    pickle.dump({'E_symbolic': E, 'P_regular': P}, open('exact_regular_hexagon.pkl', 'wb'))
