code = '''
import sys
sys.path.insert(0,'.')
sys.path.insert(0, '<project>/paper/v2/code')
import sympy as sp
import mk_integrands
from genintegrand import S_poly
t, Q = sp.symbols('t Q', positive=True)
q2 = t**2+(1-t)**2
k=5; key=(4,1)
ti = mk_integrands.t_integrand(k)
c = ti[key]
lhs = -c/(2*sp.factorial(k))
e = sp.expand(lhs)
for _ in range(12): e = sp.expand(e.replace(Q**2, q2))
qfree = sp.expand(e.subs(Q,0))
qlin = sp.expand(sp.expand((e-qfree)/Q).subs(Q,0))
print("paper coeff (M_k-normalized), Q-free part:", sp.factor(qfree))
print("paper coeff, Q-linear part (factored):", sp.factor(qlin))
m=key[0]
Sp = S_poly(k,m)
Sexpr = sum(sp.Rational(cc.numerator,cc.denominator)*t**ee for e_,cc in enumerate(Sp) for ee in [e_])
print("S_4(t) =", Sexpr)
print("claim qlin = S_4/(2t(1-t)q2^5):", sp.factor(Sexpr/(2*t*(1-t)*q2**k)))
print("ratio:", sp.simplify(qlin/(Sexpr/(2*t*(1-t)*q2**k))))
'''
open('rev_12c.py','w').write(code)
import sys, subprocess
work = '<project>/research/08_theory/WT2r_review'
r = subprocess.run([sys.executable, f'{work}/rev_12c.py'], capture_output=True, text=True, cwd=work, timeout=300)
print(r.stdout); print(r.stderr[-600:] if r.returncode else "")
