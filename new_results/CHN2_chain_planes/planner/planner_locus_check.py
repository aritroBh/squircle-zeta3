# Planner (W4) check of the 'only if' direction: for positive a, the common zeros of the weight-8 zeta(3,5) form (C3) and
# the weight-10 zeta(3,7) form (Q5, from CHN1 forms.pkl) are exactly Pi1 = {a1=a2,a3=a4} and Pi2 = {a1=a3,a2=a4}.
# (1) Q5 vanishes identically on Pi1, Pi2; (2) resultant in a4 (a3 = 1) factored; each factor's common zeros solved and
# tested for positivity; (3) random sampling of the positive part of {C3 = 0}: min |Q5| away from the planes.
import pickle, random, itertools
import sympy as sp
x, y, z, w = sp.symbols('a1 a2 a3 a4', positive=True)
X = (x, y, z, w)
forms = pickle.load(open('<project>/research/08_theory/CHN1_chain_bodies/forms.pkl', 'rb'))
def to_sp(d):
    return sum(sp.Rational(c) * sp.prod([v ** int(e) for v, e in zip(X, k.split(','))]) for k, c in d.items())
F37 = sp.expand(to_sp(forms[(10, 6)]))
F35 = sp.expand(to_sp(forms[(8, 3)]))
p, q = x + w, y + z
C3 = sp.expand(2*(p-q)*(2*p+q)*(p+2*q) - 9*x*w*(p+2*q) + 9*y*z*(2*p+q))
assert sp.expand(F35 - x*y*z*w*(x+y+z+w)*C3/10) == 0, 'P8 mismatch'
Q5 = sp.cancel(224 * F37 / (x*y*z*w*(x+y+z+w)))
Q5 = sp.expand(Q5); assert sp.Poly(Q5, *X).is_homogeneous and sp.Poly(Q5, *X).total_degree() == 5
s, t = sp.symbols('s t')
print('P8 = a1a2a3a4(Sum a)C3/10 confirmed; Q5 quintic, irreducible over Q:', len(sp.factor_list(Q5)[1]) == 1)
print('Q5 on Pi1:', sp.expand(Q5.subs({x: s, y: s, z: t, w: t})), '| on Pi2:', sp.expand(Q5.subs({x: s, y: t, z: s, w: t})))
c3 = C3.subs(z, 1); q5 = Q5.subs(z, 1)
R = sp.factor(sp.resultant(c3, q5, w))
print('Res_a4(C3,Q5)|a3=1 =', R)
for fac, mult in sp.factor_list(R)[1]:
    # solve the factor for a2 (or a1) and find common a4 roots
    sol = sp.solve(fac, y)
    if not sol:
        print(f'  factor {fac}: no a2-solution (constrains a1 only)'); continue
    for a2s in sol:
        g = sp.gcd(sp.Poly(sp.together(c3.subs(y, a2s)).as_numer_denom()[0], w), sp.Poly(sp.together(q5.subs(y, a2s)).as_numer_denom()[0], w))
        roots = sp.solve(g.as_expr(), w) if g.degree() > 0 else []
        for r in roots:
            pos = sp.solve([a2s > 0, r > 0, x > 0], x) if True else None
            print(f'  factor {fac}: a2 = {a2s}, a4 = {sp.simplify(r)}   [positive region in a1: {pos}]')
# (3) sampling
random.seed(1)
fC = sp.lambdify((x, y, z, w), C3); fQ = sp.lambdify((x, y, z, w), Q5)
best = []
for _ in range(4000):
    a1, a2, a3 = [random.uniform(0.05, 5) for _ in range(3)]
    for r in sp.Poly(C3.subs({x: a1, y: a2, z: a3}), w).nroots():
        if abs(sp.im(r)) < 1e-12 and sp.re(r) > 0:
            a4 = float(sp.re(r)); v = (a1, a2, a3, a4); nrm = sum(v)
            d = min(abs(a1 - a2) + abs(a3 - a4), abs(a1 - a3) + abs(a2 - a4)) / nrm
            best.append((d, abs(fQ(*v)) / nrm ** 5))
far = [b for b in best if b[0] > 0.05]
print(f'sampled {len(best)} positive points on C3=0; {len(far)} with plane-distance > 0.05; min |Q5|/|a|^5 among them = {min(b[1] for b in far):.3e}')
near = sorted(best)[:3]; print('closest-to-plane samples (dist, |Q5|):', [(round(d, 4), f'{qv:.2e}') for d, qv in near])
