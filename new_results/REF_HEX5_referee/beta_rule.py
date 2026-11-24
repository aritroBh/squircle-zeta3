"""Check of the Beta-function rule used for the order-five local law:
int_0^1 t^n (1-t)^k log^a t log^b(1-t) dt, a+b=5:  [zeta5] = 24 B(n+1,k+1) iff a,b>=1;
[zeta2 zeta3] = -12 B(n+1,k+1) iff (a,b) in {(3,2),(2,3)}; and no zeta5, zeta2zeta3 for a+b<=4.
Exact: sympy derivatives of Gamma(x)Gamma(y)/Gamma(x+y), polygammas at integers expanded."""
import sympy as sp
x, y = sp.symbols('x y')
Bf = sp.gamma(x) * sp.gamma(y) / sp.gamma(x + y)
z2, z3, z5 = sp.zeta(2), sp.zeta(3), sp.zeta(5)
bad = 0
for (n, k) in [(0, 0), (1, 0), (2, 3), (4, 1)]:
    for s in range(1, 6):
        for a_ in range(0, s + 1):
            b_ = s - a_
            d = sp.diff(Bf, x, a_, y, b_) if a_ and b_ else (sp.diff(Bf, x, a_) if a_ else sp.diff(Bf, y, b_))
            v = sp.expand(sp.expand_func(d.subs({x: n + 1, y: k + 1})))
            v = sp.expand(v.subs(sp.EulerGamma, 0))  # gamma cancels; check below
            Bv = sp.beta(n + 1, k + 1)
            Bv = sp.gamma(n + 1) * sp.gamma(k + 1) / sp.gamma(n + k + 2)
            c5 = v.coeff(sp.zeta(5))
            c23 = sp.expand(v).coeff(sp.zeta(3)).coeff(sp.pi, 2) * 6  # zeta2 zeta3 = pi^2 zeta3/6
            e5 = 24 * Bv if (s == 5 and a_ >= 1 and b_ >= 1) else 0
            e23 = -12 * Bv if (a_, b_) in ((3, 2), (2, 3)) else 0
            ok = sp.simplify(c5 - e5) == 0 and sp.simplify(c23 - e23) == 0
            if not ok:
                bad += 1
                print('MISMATCH', n, k, a_, b_, c5, c23)
print('beta rule mismatches:', bad)
