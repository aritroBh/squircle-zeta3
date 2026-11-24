
import sys
sys.path.insert(0,'.')
sys.path.insert(0, '<project>/paper/v2/code')
import sympy as sp
import mk_integrands
from genintegrand import S_poly
t, Q = sp.symbols('t Q')
q2 = t**2+(1-t)**2
ti4 = mk_integrands.t_integrand(4)
tv = sp.Rational(2,7)
q2v = tv**2+(1-tv)**2
for key in [(3,1),(4,0),(2,2)]:
    c = ti4[key]
    lhs = -c/(2*sp.factorial(4))
    lhs_v = lhs.subs(t, tv).subs(Q, sp.sqrt(q2v))
    lhs_v = sp.nsimplify(sp.N(lhs_v, 30))
    Sp = S_poly(4,key[0])
    Sexpr = sum(sp.Rational(cc.numerator,cc.denominator)*t**ee for ee,cc in enumerate(Sp))
    rhs_v = sp.N((sp.sqrt(q2)*Sexpr/(2*t*(1-t)*q2**4)).subs(t, tv), 30)
    print(f"key={key}: lhs={sp.N(lhs_v,12)}  rhs={sp.N(rhs_v,12)}")
    print(f"   diff = {sp.N(lhs_v-rhs_v, 5)}")
