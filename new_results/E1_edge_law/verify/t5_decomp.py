# T5: verifier's own symbolic derivation of the d^3 coefficient of |x_p'(t)| on the edges of the rational hexagon (b=0)
# directly from x_p(t) = y(t) F_p(y(t))^{-1/p}  (no use of the solver's s01 formulas), then
#   (a) the "canonical" zeta(3) part = 6*2*int(coefficient of log^2 t log(1-t) + log t log^2(1-t)) -> compare with Z_e
#   (b) exact zeta(3) content of the NON-LOCAL monomials t^n log^a t log^b(1-t) log^c(1+r't) (c>=1), via PSLQ per monomial
#   (c) consistency: local + non-local = per-arc P''' from the polar jet code
import sys, json, sympy as sp, mpmath as mp
mp.mp.dps = 60
t, d = sp.symbols('t d', positive=True)
A = sp.symbols('A0:3')   # A_j = log a_j(t)/c_j-type logs:  a_j(t) e^{d A_j} = c_j |<v_j,y>|^p
X, Y, W = sp.symbols('X Y W')   # log t, log(1-t), log(non-local affine factor normalised to 1 at t=0)
l2 = sp.log(2)

def edge_E3(V0, V1, terms):
    """terms: list of (a_j(t) [=c_j |<v_j,y(t)>| at p=1], A_j expressed via X,Y,W and logs) ; returns E3 as poly"""
    V0 = sp.Matrix(V0); V1 = sp.Matrix(V1)
    y = (1-t)*V0 + t*V1; yp = V1 - V0
    F = sum(a*sp.exp(d*Aj) for (a, Aj) in zip([q[0] for q in terms], A))
    # dF/dt: d/dt [c|<v,y>|^p] = p * a_j'(t) * e^{d A_j}   (exact, since a_j A_j' = a_j')
    Fp = (1+d)*sum(sp.diff(a, t)*sp.exp(d*Aj) for (a, Aj) in zip([q[0] for q in terms], A))
    assert sp.simplify(F.subs(d, 0) - 1) == 0
    Phi = Fp/((1+d)*F)
    G = F**(-2/(1+d))*(yp.dot(yp) - 2*y.dot(yp)*Phi + y.dot(y)*Phi**2)
    ser = sp.series(sp.sqrt(G), d, 0, 4).removeO()
    E3 = sp.expand(ser.coeff(d, 3))
    sub = {A[j]: terms[j][1] for j in range(3)}
    return sp.expand(E3.xreplace(sub))

# hexagon b=0: F_p = |y/4|^p + 3|(y-x)/4|^p + |x/4|^p   (vertices (1,0),(2,2),(0,1))
edges = {
 # edge 1: (1,0)->(2,2): y/4 = t/2, (y-x)/4 = (t-1)/4, x/4 = (1+t)/4 ; r = 2, W = log(1+t)
 1: ((1, 0), (2, 2), [(t/2, X - l2), (3*(1-t)/4, Y - 2*l2), ((1+t)/4, W - 2*l2)], 'log(1+t)'),
 # edge 3: (0,1)->(-1,0): y/4=(1-t)/4, (y-x)/4 = 1/4 (constant), x/4 -> |x|/4 = t/4 ; fully local
 3: ((0, 1), (-1, 0), [((1-t)/4, Y - 2*l2), (3*sp.Rational(1, 4), -2*l2 + 0*W), (t/4, X - 2*l2)], None),
}
from functools import lru_cache
@lru_cache(None)
def beta_moment(n, a, b):
    x, yv = sp.symbols('x yv')
    Bf = sp.gamma(x)*sp.gamma(yv)/sp.gamma(x+yv)
    expr = sp.diff(Bf, x, a, yv, b) if (a or b) else Bf
    val = expr.subs({x: n+1, yv: 1})
    val = sp.expand(sp.expand_func(sp.simplify(val)))
    return val

basis = lambda: [mp.zeta(3), mp.pi**2*mp.log(2), mp.log(2)**3, mp.pi**2, mp.log(2)**2, mp.log(2), mp.mpf(1)]
res = {}
for key, (V0, V1, terms, Wname) in edges.items():
    E3 = edge_E3(V0, V1, terms)
    P = sp.Poly(E3, t, X, Y, W)
    loc = 0; nonloc_terms = []
    for (n, a, b, c), co in P.terms():
        if c == 0:
            loc += co*beta_moment(n, a, b)
        else:
            nonloc_terms.append((n, a, b, c, co))
    loc = sp.expand(6*loc)
    z3_loc = sp.nsimplify(sp.simplify(loc.coeff(sp.zeta(3))))
    print(f"== edge {key}: 6*int(local part) zeta(3) coefficient = {z3_loc}  (= {sp.N(z3_loc, 25)})")
    print(f"   local part (exact) = {sp.simplify(loc)}")
    # non-local numerics + PSLQ identification monomial by monomial
    z3_nl = 0; nl_val = mp.mpf(0)
    for (n, a, b, c, co) in nonloc_terms:
        f = lambda tt: tt**n*mp.log(tt)**a*mp.log(1-tt)**b*mp.log(1+tt)**c
        I = mp.quad(f, [0, 0.5, 1])
        rel = mp.pslq([I] + basis(), maxcoeff=10**8, maxsteps=10**6)
        assert rel and rel[0] != 0, (n, a, b, c, rel)
        cz3 = sp.Rational(-rel[1], rel[0])
        z3_nl += co*cz3
        nl_val += mp.mpf(sp.N(co, 70))*I
    z3_nl = sp.nsimplify(sp.simplify(6*z3_nl))
    print(f"   non-local monomials: {len(nonloc_terms)}; their exact zeta(3) content (x6) = {z3_nl} (= {sp.N(z3_nl, 25)})")
    tot = mp.mpf(sp.N(loc.subs(sp.zeta(3), sp.zeta(3)), 70)) + 6*nl_val
    print(f"   6*int E3 (local exact + non-local quad) = {mp.nstr(tot, 40)}")
    print(f"   total zeta(3) coefficient = {sp.nsimplify(z3_loc + z3_nl)} (= {sp.N(z3_loc + z3_nl, 25)})")
    res[key] = dict(z3_local=str(z3_loc), z3_nonlocal=str(z3_nl), total=mp.nstr(tot, 50), local=str(loc))
json.dump(res, open('research/10_frontier/E1_edge_law/verify/out/t5_decomp.json', 'w'), indent=1)
