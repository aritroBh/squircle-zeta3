WD='<project>/research/08_theory/FLD1_field_rule_proof'
with open(WD+'/progress.log','a') as fh:
    fh.write("step23: Z-sum canonical form mu_zsum_K5.pkl; k=3 reproduces mu_zsum.pkl (B3 template) EXACTLY\n")
code = r'''"""step24_inspect.py -- inspect the Z-sum coefficients rho_w(n) of mu_{k,n}: proper? poles? polynomial part?"""
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
'''
open(WD+'/step_24_inspect.py','w').write(code)
import subprocess
VENV='<project>/.venv/bin/python'
r = subprocess.run(['nice',VENV,'step_24_inspect.py'],capture_output=True,text=True,cwd=WD)
print(r.stdout); print(r.stderr[-2000:] if r.returncode else '')