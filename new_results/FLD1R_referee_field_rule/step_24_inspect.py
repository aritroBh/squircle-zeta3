"""step24_inspect.py -- inspect the Z-sum coefficients rho_w(n) of mu_{k,n}: proper? poles? polynomial part?"""
import sys, pickle
import sympy as sp
n = sp.Symbol('n')
Z = pickle.load(open('mu_zsum_K5.pkl','rb'))
def analyze(k):
    print(f'--- k={k} ---')
    maxdeg_poly = -999
    poles = set()
    for zk, terms in Z[k].items():
        for w, cs in terms:
            e = sp.together(sp.sympify(cs))
            num, den = sp.fraction(e)
            dn = sp.Poly(sp.expand(num), n)
            dd = sp.Poly(sp.expand(den), n)
            degpoly = dn.degree() - dd.degree()
            maxdeg_poly = max(maxdeg_poly, degpoly)
            # roots of den
            for root in sp.roots(sp.expand(den), n):
                poles.add(sp.nsimplify(root))
    print(f'   max (deg num - deg den) = {maxdeg_poly}   (<=0 means proper, no polynomial part)')
    print(f'   poles of coefficients: {sorted(poles, key=sp.N)}')
for k in (1,2,3,4,5):
    analyze(k)
