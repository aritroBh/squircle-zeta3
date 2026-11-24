
ctrl = r'''
# ctrl.py : area line (p=r) exact coefficients as Gamma-quantity control
import sympy as sp, pickle
s1 = sp.symbols('s1')
p = 2/(1+2*s1)
area = 4*sp.gamma(1+1/p)**2/sp.gamma(1+2/p)
# expand via loggamma series manually
lg = 2*sp.loggamma(1+1/p) - sp.loggamma(1+2/p) + sp.log(4)
lgser = sp.series(lg, s1, 0, 5).removeO()
lgser = sp.expand_func(lgser)
lgser = sp.expand(lgser)
aser = sp.exp(sp.series(lgser, s1, 0, 5).removeO())
aser = sp.series(aser, s1, 0, 5).removeO().expand()
print("area(K_p) expansion at p=2 (control: a Gamma quantity):")
for n in range(0,5):
    print(f"  s^{n}:", sp.simplify(aser.coeff(s1,n)), "  =", sp.N(aser.coeff(s1,n),20))
# compare with bivariate area-line sums
b = {eval(k): sp.sympify(v) for k,v in pickle.load(open('scratchpad/circle_bivariate.pkl','rb')).items()}
print("\nbivariate area-line sums vs exact area coefficients:")
for n in range(1,5):
    val = sp.expand(sum(b[(j,n-j)] for j in range(n+1)))
    ref = sp.expand_func(sp.expand(aser.coeff(s1,n)))
    ref = sp.expand(ref.subs(sp.polygamma(1, sp.Rational(3,2)), sp.pi**2/2-4).subs(sp.polygamma(2, sp.Rational(3,2)), 16-14*sp.zeta(3)))
    print(f"  order {n}: diff = {sp.simplify(sp.expand(val-ref))}")
'''
open('ctrl.py','w').write(ctrl)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','ctrl.py'],capture_output=True,text=True,timeout=900)
print(r.stdout); print(r.stderr[-400:])
