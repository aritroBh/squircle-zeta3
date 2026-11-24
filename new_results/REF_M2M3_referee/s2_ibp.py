"""Referee step 2a: exact re-derivation (own code, sympy) of the weight split
   M_2 = W_3 + T,   M_3 = W_4 + W_3' + T
from the definition M_k = -(1/2)[eps^k] int_0^1 Phi dt, via the Euler substitution,
folding v -> 1/v, and integration by parts of the poles at v = +-i.

Kernels are tracked as a dict {(i,j): R(v)} meaning int_b^a R(v) lt^i lmt^j dv.
"""
import sympy as sp

v, e = sp.symbols('v e')
x, y = sp.symbols('x y')          # x = log t, y = log(1-t)
A, B, S = sp.symbols('A B S', positive=True)
s2 = sp.sqrt(2)
a = 1 + s2
b = s2 - 1

# ---- 1. eps-expansion of Phi = A e^{e x} + B e^{e y} - sqrt(A^2 e^{2ex} + B^2 e^{2ey}),  S = sqrt(A^2+B^2)
u = (A**2*(sp.exp(2*e*x) - 1) + B**2*(sp.exp(2*e*y) - 1))/S**2
Phi = A*sp.exp(e*x) + B*sp.exp(e*y) - S*sp.sqrt(1 + u)
ser = sp.expand(sp.series(Phi, e, 0, 4).removeO())

tv = (v - b)*(v + a)/(4*v)
qv = (v**2 + 1)/(2*s2*v)
dtdv = (v**2 + 1)/(4*v**2)
subsd = {A: 1/tv, B: 1/(1 - tv), S: qv/(tv*(1 - tv))}


def canon(R):
    return sp.factor(sp.cancel(sp.together(sp.radsimp(R))), extension=s2)


def kernels(k):
    ck = sp.expand(ser.coeff(e, k))
    P = sp.Poly(ck, x, y)
    K = {}
    for (i, j), c in zip(P.monoms(), P.coeffs()):
        # M_k = -1/2 int Phi_k dt
        R = canon(-sp.Rational(1, 2)*c.subs(subsd)*dtdv)
        K[(i, j)] = sp.simplify(K.get((i, j), 0) + R)
    return K


def fold(K):
    out = {}
    for (i, j), R in K.items():
        if i >= j:
            out[(i, j)] = canon(out.get((i, j), 0) + R)
        else:
            out[(j, i)] = canon(out.get((j, i), 0) + R.subs(v, 1/v)/v**2)
    return {kk: R for kk, R in out.items() if sp.simplify(R) != 0}


dlt = (v**2 + 1)/(v*(v - b)*(v + a))        # d log t / dv
dlmt = -(v**2 + 1)/(v*(a - v)*(v + b))      # d log(1-t) / dv


def split_ipoles(R):
    """R = Rreal + Rimag where Rimag collects the partial fractions with (v^2+1)^n denominators."""
    num, den = sp.fraction(sp.together(R))
    terms = sp.Add.make_args(sp.apart(sp.together(R), v, extension=s2))
    Ri = 0
    for tterm in terms:
        _, d = sp.fraction(sp.together(tterm))
        if sp.Poly(d, v).rem(sp.Poly(v**2 + 1, v)).is_zero:
            Ri += tterm
    return sp.together(R - Ri), sp.together(Ri)


def herm(Ri):
    """Ri = sum (c v + d)/(v^2+1)^n. Return rational G and simple-pole rest with Ri = G' + rest."""
    G = 0
    rest = 0
    # coefficients by power n
    coeffs = {}
    for tt in sp.Add.make_args(sp.expand(Ri)):
        num, den = sp.fraction(sp.together(tt))
        n = 0
        dp = sp.Poly(den, v)
        q = sp.Poly(v**2 + 1, v)
        cst = dp
        while cst.rem(q).is_zero:
            cst = cst.quo(q)
            n += 1
        assert cst.degree() == 0
        pn = sp.Poly(sp.expand(num/cst.as_expr()), v)
        # reduce numerator modulo (v^2+1) repeatedly
        while pn.degree() >= 2:
            qq, rr = pn.div(q)
            coeffs[n] = coeffs.get(n, [0, 0])
            c1 = rr.coeff_monomial(v); c0 = rr.coeff_monomial(1)
            coeffs[n][0] += c1; coeffs[n][1] += c0
            pn = qq; n -= 1
        coeffs[n] = coeffs.get(n, [0, 0])
        coeffs[n][0] += pn.coeff_monomial(v); coeffs[n][1] += pn.coeff_monomial(1)
    assert all(n >= 1 for n in coeffs)
    for n in sorted(coeffs, reverse=True):
        c, d = coeffs[n]
        if n == 1:
            rest += (c*v + d)/(v**2 + 1)
            continue
        # int v/(v^2+1)^n = -1/(2(n-1)(v^2+1)^(n-1))
        G += -c/(2*(n - 1)*(v**2 + 1)**(n - 1))
        # int 1/(v^2+1)^n = v/(2(n-1)(v^2+1)^(n-1)) + (2n-3)/(2(n-1)) int 1/(v^2+1)^(n-1)
        G += d*v/(2*(n - 1)*(v**2 + 1)**(n - 1))
        coeffs[n - 1] = coeffs.get(n - 1, [0, 0])
        coeffs[n - 1][1] += d*sp.Rational(2*n - 3, 2*(n - 1))
    assert sp.simplify(sp.diff(G, v) + rest - Ri) == 0
    return sp.together(G), rest


def ibp(K):
    """Remove poles at +-i of order >= 2 by parts, from top degree down. Returns final kernels and
    the residual simple-pole (+-i) parts, which must be zero."""
    K = fold(K)
    simple = {}
    maxdeg = max(i + j for i, j in K)
    for d in range(maxdeg, 0, -1):
        K = fold(K)
        for (i, j) in [kk for kk in K if sum(kk) == d]:
            Rr, Ri = split_ipoles(K[(i, j)])
            if Ri == 0:
                continue
            G, rem = herm(Ri)
            if sp.simplify(rem) != 0:
                simple[(i, j)] = sp.simplify(rem)
            G = sp.together(G - G.subs(v, b))   # normalisation G(b) = 0
            assert sp.simplify(sp.radsimp(G.subs(v, b))) == 0
            K[(i, j)] = canon(Rr + (Ri - sp.diff(G, v)))
            # int G' L = [G L]_b^a - int G dL ; boundary terms vanish (i >= 1; G(b) = 0)
            if i >= 1:
                K[(i - 1, j)] = canon(K.get((i - 1, j), 0) - G*i*dlt)
            if j >= 1:
                K[(i, j - 1)] = canon(K.get((i, j - 1), 0) - G*j*dlmt)
    return fold(K), simple


def sym_diag(K):
    out = {}
    for (i, j), R in K.items():
        if i == j:
            out[(i, j)] = canon((R + R.subs(v, 1/v)/v**2)/2)
        else:
            out[(i, j)] = R
    return out


printed = {
    2: {(2, 0): (4 - 3*s2)/(8*v) - 1/(v + a), (1, 1): s2/8/v,
        (1, 0): 1/(2*v) + 1/(2*(v - a)) - 1/(v + a) - 1/(2*(v + b))},
    3: {(3, 0): (8 - 7*s2)/(48*v) - 1/(3*(v + a)), (2, 1): s2/16/v,
        (2, 0): (1 - s2)/(4*v) + 1/(4*(v - a)) - 1/(4*(v + a)) - 1/(2*(v + b)),
        (1, 1): (1 + s2)/(4*v) - 1/(2*(v + a)),
        (1, 0): 1/(2*v) + 1/(2*(v - a)) - 1/(v + a) - 1/(2*(v + b))},
}

if __name__ == '__main__':
    for k in (1, 2, 3):
        K0 = kernels(k)
        print(f'=== k={k}: raw kernels after Euler substitution (before folding)')
        for kk, R in sorted(K0.items(), reverse=True):
            print('  ', kk, R)
        Kf, simple = ibp(K0)
        Kf = sym_diag(Kf)
        print(f'--- k={k}: after fold + IBP (diagonal symmetrised)')
        for kk, R in sorted(Kf.items(), reverse=True):
            print('  ', kk, sp.apart(R, v, extension=s2))
        print('   residual simple poles at +-i:', simple if simple else 'none')
        if k in printed:
            P = sym_diag({kk: canon(R) for kk, R in printed[k].items()})
            keys = set(P) | set(Kf)
            for kk in sorted(keys, reverse=True):
                diff = canon(Kf.get(kk, 0) - P.get(kk, 0))
                print('   kernel', kk, 'own - printed =', sp.simplify(diff))
