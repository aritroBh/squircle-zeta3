# g9_local5.py -- exact P^(5)(1) of fully local polygons, from the length element of Lemma g4:lem:arc (own code).
# Edge data: terms a_j(t) = alpha t, beta (1-t), and constant far terms rho_k with logarithms lam_k;
# G = F^{-1/p} sqrt((L - sigma Phi)^2 + h^2 Phi^2),  F = sum a_j e^{d lam_j},  D = sum a_j' e^{d lam_j},  Phi = D/F,  p = 1+d.
# Series in d to order 5 with coefficients polynomial in t, X = log t, Y = log(1-t); integrals by the formal Beta moments of g9_beta5.py.
# Examples: diamond (m=2), regular hexagon (closed form of the paper, 15 terms), sheared parallelogram +-(2,0), +-(1,2).
import sympy as sp, mpmath as mp, time, sys
sys.path.insert(0, '.')
import g9_beta5 as bl
from sympy import Rational as R
t, X, Y, eta = sp.symbols('t X Y eta')
z2, z3, z5 = bl.z2, bl.z3, bl.z5
K = 6
def mul(a, b): return [sp.expand(sum(a[i]*b[k-i] for i in range(k+1))) for k in range(K)]
def inv(a):
    r = [sp.Integer(1)] + [0]*(K-1)
    for k in range(1, K): r[k] = sp.expand(-sum(a[i]*r[k-i] for i in range(1, k+1)))
    return r
def logs(a):          # log of series with a[0]=1
    x = [0] + a[1:]; res = [0]*K; pw = [sp.Integer(1)] + [0]*(K-1)
    for j in range(1, K):
        pw = mul(pw, x); res = [sp.expand(res[k] + R((-1)**(j+1), j)*pw[k]) for k in range(K)]
    return res
def expo(a):          # exp of series with a[0]=0
    res = [sp.Integer(1)] + [0]*(K-1); pw = [sp.Integer(1)] + [0]*(K-1)
    for j in range(1, K):
        pw = mul(pw, a); res = [sp.expand(res[k] + pw[k]/sp.factorial(j)) for k in range(K)]
    return res
def edge_series(alpha, beta, far, mu0, mu1, u0, m):
    """returns [E_0..E_5]/L as polynomials in t,X,Y (and symbols in mu's / far logs)."""
    lam = [(alpha*t, X+mu0), (beta*(1-t), Y+mu1)] + [(rho, l) for rho, l in far]
    dlam = [(alpha, X+mu0), (-beta, Y+mu1)]
    F = [sp.Integer(1)] + [sp.expand(sum(a*l**k for a, l in lam)/sp.factorial(k)) for k in range(1, K)]
    D = [sp.expand(sum(a*l**k for a, l in dlam)/sp.factorial(k)) for k in range(K)]
    Phi = mul(D, inv(F))
    Fpow = expo([sp.expand(-c) for c in mul([(-1)**k for k in range(K)], logs(F))])    # F^{-1/(1+d)}
    Xi = [sp.expand(-(u0+t)*c) for c in Phi]; Xi[0] += 1                                  # 1 - (sigma/L) Phi
    Q = mul(Xi, Xi); Phi2 = mul(Phi, Phi)
    Q = [sp.expand(Q[k] + m*Phi2[k]) for k in range(K)]                                    # (..)^2 + m Phi^2
    x = [0] + Q[1:]; sq = [sp.Integer(1)] + [0]*(K-1); pw = [sp.Integer(1)] + [0]*(K-1)
    for j in range(1, K):
        pw = mul(pw, x); sq = [sp.expand(sq[k] + sp.binomial(R(1, 2), j)*pw[k]) for k in range(K)]
    return mul(Fpow, sq)
_cache = {}
def moment(n, a, b):
    if (n,) not in _cache: _cache[(n,)] = bl.Ihat(n, 0)[0]
    return _cache[(n,)][(a, b)]
def integrate(poly):
    poly = sp.Poly(sp.expand(poly), t, X, Y)
    tot = 0
    for (n, a, b), c in poly.terms(): tot += c*moment(n, a, b)
    return sp.expand(tot)
def edge_integrals(*args):
    E = edge_series(*args)
    return [integrate(e) for e in E]       # int E_k dt / L
t0 = time.time()
# ---- diamond: 4 edges, L = sqrt2, u0 = -1/2, m = 1/4, alpha = beta = 1
I = edge_integrals(1, 1, [], 0, 0, -R(1, 2), R(1, 4))
P5 = 120*4*sp.sqrt(2)*I[5]
want = sp.sqrt(2)*(-90*z5 - 210*z2*z3 + 396*z2**2 + 720*z3 - 480*z2)
print('diamond: P5(1) - sqrt2(-90 z5 - 210 z2 z3 + 396 z2^2 + 720 z3 - 480 z2) =', sp.expand(P5 - want), '  (%.0f s)' % (time.time()-t0))
print('         lower orders: P(1)=%s, P\'(1)=%s' % (sp.expand(4*sp.sqrt(2)*I[0]), sp.expand(4*sp.sqrt(2)*I[1])))
# ---- regular hexagon: 6 edges, L = 1/sqrt3, u0 = -1/2, m = 3/4, alpha = beta = rho = 1/2, logs -eta
E5mon = edge_series(R(1, 2), R(1, 2), [(R(1, 2), -eta)], -eta, -eta, -R(1, 2), R(3, 4))
print('regular hexagon: number of monomials t^n log^a t log^b(1-t) in the d^5 integrand of an edge:', len(sp.Poly(E5mon[5], t, X, Y).terms()), ' (the paper says 91)')
I = edge_integrals(R(1, 2), R(1, 2), [(R(1, 2), -eta)], -eta, -eta, -R(1, 2), R(3, 4))
P5 = 120*6*I[5]/sp.sqrt(3)
want = sp.sqrt(3)/128*(30720*eta - 4920*6*z2 - 61440*eta**2 + 6240*6*z2*eta + 30720*eta**3 + 448*36*z2**2
      + 5760*z3*eta - 1920*6*z2*eta**2 - 5120*eta**4 + 11340*z5 + 5580*z2*z3 - 123*36*z2**2*eta
      - 1440*z3*eta**2 + 160*6*z2*eta**3 + 256*eta**5)
print('regular hexagon: P5(1) - (paper closed form) =', sp.expand(P5 - want), '  (%.0f s)' % (time.time()-t0))
mp.mp.dps = 40
val = sp.N(want.subs({z2: sp.pi**2/6, z3: sp.zeta(3), z5: sp.zeta(5), eta: sp.log(2)}), 40)
print('         numerical value of the closed form:', val)
print('         [z5] =', sp.simplify(sp.Poly(sp.expand(P5), z2, z3, z5, eta).coeff_monomial(z5)), ' [z2 z3] =', sp.simplify(sp.Poly(sp.expand(P5), z2, z3, z5, eta).coeff_monomial(z2*z3)))
# ---- sheared parallelogram +-(2,0), +-(1,2): two edges of length sqrt5 (u0=-2/5, m=16/25) and two of length sqrt13 (u0=-7/13, m=16/169)
Ia = edge_integrals(1, 1, [], 0, 0, -R(2, 5), R(16, 25))
Ib = edge_integrals(1, 1, [], 0, 0, -R(7, 13), R(16, 169))
P5 = 120*2*(sp.sqrt(5)*Ia[5] + sp.sqrt(13)*Ib[5])
s5, s13 = sp.sqrt(5), sp.sqrt(13)
want = (z5*(-6912*s5/25 - 172800*s13/28561) + z2*z3*(-10752*s5/25 - 821760*s13/28561) + z2**2*(119808*s5/125 + 1262592*s13/28561)
        + z3*(4608*s5/5 + 23040*s13/169) + z2*(-3072*s5/5 - 15360*s13/169))
print('sheared parallelogram: P5(1) - (value in our notes) =', sp.expand(P5 - want), '  (%.0f s)' % (time.time()-t0))
