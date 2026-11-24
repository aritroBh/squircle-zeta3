code = '''
import sys
sys.path.insert(0,'.')
sys.path.insert(0, '<project>/paper/v2/code')
import sympy as sp
import mk_integrands
from genintegrand import S_poly
t, Q = sp.symbols('t Q')
q2 = t**2+(1-t)**2
ti4 = mk_integrands.t_integrand(4)
c = ti4[(3,1)]
print("raw c:", sp.factor(c))
lhs = sp.together(-c/(2*sp.factorial(4)))
N, D = sp.fraction(lhs)
print("N:", sp.factor(N)); print("D:", sp.factor(D))
# evaluate numerically with Q -> sqrt(q2) at t=2/7, exact via nsimplify on algebraic
tv = sp.Rational(2,7)
q2v = tv**2+(1-tv)**2
import mpmath
qval = mpmath.sqrt(float(q2v))
lhs_v = float(lhs.subs({t: tv, Q: qval}))
Sp = S_poly(4,3)
Sexpr = sum(sp.Rational(cc.numerator,cc.denominator)*t**ee for ee,cc in enumerate(Sp))
rhs_v = float((sp.sqrt(q2)*Sexpr/(2*t*(1-t)*q2**4)).subs(t, tv))
print("lhs(=-c/48) at t=2/7:", lhs_v)
print("rhs(=q S_3/(2t(1-t)q2^4)) at t=2/7:", rhs_v)
print("ratio:", lhs_v/rhs_v)
# also try with an extra minus
print("rhs with - sign:", -rhs_v)
# and maybe I mislabeled: try S_1(t) (the l2^3 l1 piece?)
S1 = S_poly(4,1)
S1e = sum(sp.Rational(cc.numerator,cc.denominator)*t**ee for ee,cc in enumerate(S1))
rhs1_v = float((sp.sqrt(q2)*S1e/(2*t*(1-t)*q2**4)).subs(t, tv))
print("rhs with S_1:", rhs1_v)
'''
open('rev_12g.py','w').write(code)
import sys, subprocess
work = '<project>/research/08_theory/WT2r_review'
r = subprocess.run([sys.executable, f'{work}/rev_12g.py'], capture_output=True, text=True, cwd=work, timeout=300)
print(r.stdout); print(r.stderr[-600:] if r.returncode else "")
