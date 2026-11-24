"""Referee's own exact derivation of the order-five local law (prop:order5).

Edge y(t) = V0 + t D, |D|^2 = L^2, V0.D = q0, |V0|^2 = V02.  Terms on the edge at p = 1:
A0 t (line through V0), A1 (1-t) (line through V1), and the rest, summing to 1.
At p = 1+d the j-th term is w_j^{-d} l_j^{1+d}.  zeta(5) and zeta(2)zeta(3) of int t^n log^a t log^b(1-t)
need a + b = 5 with a, b >= 1 (Beta-function rule, re-derived in beta_rule.py), and at order d^5 such a
monomial is "leading-log": each d carries log t or log(1-t).  So its coefficient is obtained by replacing
t^d -> E0 = e^{s0}, (1-t)^d -> E1 = e^{s1}, setting every other d to 0, and taking [s0^a s1^b] of the arc
length element |d/dt (y / F)|, with F = 1 + A0 t (E0-1) + A1 (1-t)(E1-1) and F' = A0 (E0-1) - A1 (E1-1)
(the remaining terms add up to 1 - A0 t - A1 (1-t) and carry no leading log).
Then [zeta5] = 5! * 24 * sum_{a+b=5, a,b>=1} int_0^1 g_ab,  [zeta2 zeta3] = -5! * 12 * int (g_32 + g_23).
"""
import sympy as sp

t, eps, a, b = sp.symbols('t epsilon a b')
A0, A1, L, L0, h = sp.symbols('A0 A1 L L0 h', positive=True)


def law():
    from exact_local_edge import smul, sexp_nc, slog1p, ssqrt1p, NO
    q0 = -L0 * L
    L2 = L ** 2
    V02 = h ** 2 + L0 ** 2
    # series in eps (orders 0..5); E0 - 1 = e^{eps a} - 1, E1 - 1 = e^{eps b} - 1
    E0m = [sp.Integer(0)] + [a ** k / sp.factorial(k) for k in range(1, NO)]
    E1m = [sp.Integer(0)] + [b ** k / sp.factorial(k) for k in range(1, NO)]
    u = [sp.expand(A0 * t * E0m[k] + A1 * (1 - t) * E1m[k]) for k in range(NO)]      # F - 1
    Fp = [sp.expand(A0 * E0m[k] - A1 * E1m[k]) for k in range(NO)]
    invF = sexp_nc([-c for c in slog1p(u)])
    g = smul(Fp, invF)                                                                  # F'/F
    yD = q0 + t * L2
    yy = V02 + 2 * t * q0 + t ** 2 * L2
    gg = smul(g, g)
    # |x'|^2 = (L2 - 2 yD g + yy g^2) / F^2  ->  |x'| = L * sqrt(1 + inner) / F
    inner = [sp.expand((-2 * yD * g[k] + yy * gg[k]) / L2) for k in range(NO)]
    arc = smul(ssqrt1p(inner), invF)
    c5 = sp.expand(L * arc[5])
    P = sp.Poly(c5, a, b)
    g5 = {}
    for i in range(6):
        g5[(i, 5 - i)] = sp.expand(P.coeff_monomial(a ** i * b ** (5 - i)))
    I = {k: sp.integrate(v, (t, 0, 1)) for k, v in g5.items()}
    Z5 = sp.factor(120 * 24 * sum(I[(i, 5 - i)] for i in range(1, 5)))
    Z23 = sp.factor(-120 * 12 * (I[(3, 2)] + I[(2, 3)]))
    return Z5, Z23, g5


def printed():
    L1 = L - L0
    xb = (A0 * L0 + A1 * L1) / L
    xp = A0 + A1
    xx = A0 * A1
    Z5 = 360 * xx * h ** 2 / L * (20 * xb ** 3 - 32 * xb ** 2 + 15 * xb - 2 - 4 * (5 * xb - 4) * xx * L0 * L1 / L ** 2
                                 + h ** 2 / L ** 2 * (5 * xx * xp + 8 * xp ** 2 - 28 * xx - 15 * (xp ** 2 - 3 * xx) * xb))
    Z23 = 120 * xx * h ** 2 / L * (12 * xb ** 2 - 12 * xb + 2 - 12 * (5 * xb - 4) * xx * L0 * L1 / L ** 2
                                  + h ** 2 / L ** 2 * (15 * xx * xp - 3 * xp ** 2 + 24 * xx - 45 * xx * xb))
    return Z5, Z23


if __name__ == '__main__':
    import time
    t0 = time.time()
    Z5, Z23, g = law()
    P5, P23 = printed()
    print('derived Z5   =', Z5)
    print('derived Z23  =', Z23)
    print('Z5  derived - printed =', sp.simplify(sp.expand(Z5 - P5)))
    print('Z23 derived - printed =', sp.simplify(sp.expand(Z23 - P23)))
    # pure-log terms (a or b = 0) give no zeta5; check g_50, g_05 exist (they integrate to rationals)
    print('time %.1f s' % (time.time() - t0))
    # parallelogram (A0 = A1 = 1): printed short forms
    u_ = L0 / L
    par5 = 360 * h ** 2 / L * (1 - 4 * u_ * (1 - u_) - h ** 2 / L ** 2)
    par23 = 120 * h ** 2 / L * (2 - 12 * u_ * (1 - u_) - 3 * h ** 2 / L ** 2)
    print('parallelogram Z5  - printed short form:', sp.simplify(Z5.subs({A0: 1, A1: 1}) - par5))
    print('parallelogram Z23 - printed short form:', sp.simplify(Z23.subs({A0: 1, A1: 1}) - par23))
    import pickle
    pickle.dump({'Z5': Z5, 'Z23': Z23}, open('local_law_derived.pkl', 'wb'))
