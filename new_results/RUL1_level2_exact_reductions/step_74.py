import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import pickle, sympy as sp
WD = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
coef_acc = pickle.load(open(f'{WD}/coef_acc_raw.pkl','rb'))
print("monomials:", len(coef_acc))
items = list(coef_acc.items())
for (at,xe,ze),co in items[:6]:
    print("at", at, "xe", xe, "ze", ze, "co", co)
A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
l2 = sp.Symbol('l2')
rd1 = {A0:sp.Rational(1,2),A1:sp.Rational(3,4),B:sp.Rational(1,4),kap:sp.Integer(1),
      L2:sp.Integer(5),q0:sp.Integer(1),V02:sp.Integer(1),lA0:-l2,lA1:-2*l2,lB:-2*l2}
# sample coefficient evaluation and total log-power distribution
maxpow = 0
for (at,xe,ze),co in items:
    cov = sp.expand(co.subs(rd1))
    if cov.has(sp.sqrt(5)):
        pass
    p = sp.Poly(cov, l2) if cov != 0 else None
    if p: maxpow = max(maxpow, p.degree())
print("max l2-power in coefficients:", maxpow)
# check types of coefficients: any sqrt?
s = set()
for (at,xe,ze),co in items:
    for a in sp.preorder_traversal(co):
        s.add(type(a).__name__)
print(s)
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=600)
print(r.stdout); print(r.stderr[-600:])
