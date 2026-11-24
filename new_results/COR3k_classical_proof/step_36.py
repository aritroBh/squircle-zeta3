add = '''

# ---- proved Li2 atom reductions (classical: Euler/Landen/five-term instances) ----
LI2_REDUCE = {
    'x': sp.polylog(2, x),
    '-x': sp.polylog(2, -x),
    '1-x': -sp.polylog(2, x) + P2/6 - logx*l1m,
    '1/(1+x)': sp.polylog(2, -x) + P2/6 + logx*l1p - l1p**2/2,
    '(1+x)/2': sp.polylog(2, (1+x)/2),
    '(1-x)/(1+x)': sp.polylog(2, (1-x)/(1+x)),
    '-(1-x)/(1+x)': sp.polylog(2, (1+x)/2) - P2/6 + l1m*l1p - l1m*L2 - l1p**2/2 + L2**2/2,
    '2x/(1+x)': -sp.polylog(2, (1-x)/(1+x)) + P2/6 - (L2 + logx - l1p)*(l1m - l1p),
}
LI2_BASIS = [sp.polylog(2, x), sp.polylog(2, -x), sp.polylog(2, (1 + x)/2),
             sp.polylog(2, (1 - x)/(1 + x))]

def reduce_li2(expr):
    """replace every polylog(2, a) (a among the 8 args) by its PROVED basis form."""
    rep = {}
    for pl in sp.preorder_traversal(expr):
        if isinstance(pl, sp.polylog) and pl.args[0] == 2:
            name = arg_key(pl.args[1])
            assert name is not None, f"Li2 with unknown argument {pl.args[1]}"
            rep[pl] = LI2_REDUCE[name]
    return expr.subs(rep, simultaneous=True)

def check_weight2_zero(expr):
    """expr (after reduce_li2) must be a sum  rational(x) * (monomial in the 4
    Li2 basis atoms and the 4 log atoms).  Check every coefficient is 0."""
    e = sp.expand(reduce_li2(expr))
    gens = LI2_BASIS + list(_YS)
    es = e.subs(dict(zip(LOG_ATOMS, _YS)))
    p = sp.Poly(es, *gens, domain=sp.QQ.frac_field(x))
    terms = []
    ok = True
    for monom, coeff in p.terms():
        c = sp.cancel(sp.sympify(coeff))
        terms.append((monom, c))
        if c != 0:
            ok = False
    return ok, terms
'''
src = open("cor3_common.py").read()
assert "def check_weight1_zero" in src and "reduce_li2" not in src
open("cor3_common.py", "w").write(src + add)
print("added reduce_li2 / check_weight2_zero")
