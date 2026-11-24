# Planner (W4): exact component analysis of {C3 = 0} ∩ {Q5 = 0} (a3 = 1), no sign assumptions on symbols.
# For each linear factor of Res_a4(C3,Q5): parametrise it, take gcd in a4 of C3 and Q5 over Q(param), list the common
# components (a1,a2,1,a4) and decide whether they meet the open positive orthant.
import pickle
import sympy as sp
x, y, z, w, t = sp.symbols('a1 a2 a3 a4 t')
X = (x, y, z, w)
forms = pickle.load(open('<project>/research/08_theory/CHN1_chain_bodies/forms.pkl', 'rb'))
to_sp = lambda d: sum(sp.Rational(c) * sp.prod([v ** int(e) for v, e in zip(X, k.split(','))]) for k, c in d.items())
p, q = x + w, y + z
C3 = sp.expand(2*(p-q)*(2*p+q)*(p+2*q) - 9*x*w*(p+2*q) + 9*y*z*(2*p+q))
Q5 = sp.expand(sp.cancel(224 * sp.expand(to_sp(forms[(10, 6)])) / (x*y*z*w*(x+y+z+w))))
c3, q5 = C3.subs(z, 1), Q5.subs(z, 1)
R = sp.resultant(c3, q5, w)
for fac, _ in sp.factor_list(R)[1]:
    # parametrise the factor's zero set in the (a1,a2) plane by t
    if fac.has(y):
        a1v, a2v = t, sp.solve(fac, y)[0].subs(x, t)
    else:
        a1v, a2v = sp.solve(fac, x)[0], t
    g = sp.gcd(sp.Poly(sp.expand(c3.subs({x: a1v, y: a2v})), w), sp.Poly(sp.expand(q5.subs({x: a1v, y: a2v})), w))
    comps = sp.solve(g.as_expr(), w) if g.degree() > 0 else []
    for a4v in comps:
        pt = (sp.simplify(a1v), sp.simplify(a2v), 1, sp.simplify(a4v))
        # positivity: is there t with all four coordinates > 0 ?
        conds = [sp.simplify(c) for c in pt if sp.simplify(c).has(t)]
        const = [c for c in pt if not sp.sympify(c).has(t)]
        ok_const = all(sp.sympify(c) > 0 for c in const)
        reg = sp.reduce_inequalities([c > 0 for c in conds], t) if conds else sp.true
        tag = ''
        if pt[0] == pt[1] and pt[2] == pt[3]: tag = '  <- Pi1'
        if pt[0] == pt[2] and pt[1] == pt[3]: tag = '  <- Pi2'
        print(f'factor {fac}: component (a1,a2,a3,a4) = {pt}; positive for t in: {reg if ok_const else False}{tag}')
