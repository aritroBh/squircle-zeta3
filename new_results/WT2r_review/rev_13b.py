
import sys
sys.path.insert(0,'.')
sys.path.insert(0, '<project>/paper/v2/code')
import sympy as sp
import mk_integrands
from genintegrand import S_poly
t = mk_integrands.t; Qs = mk_integrands.Qs
q2 = t**2+(1-t)**2
k=4
ti = mk_integrands.t_integrand(k)
c = ti[(3,1)]
lhs = -c/(2*sp.factorial(k))
tv = sp.Rational(2,7); q2v = tv**2+(1-tv)**2
lhs_v = sp.N(lhs.subs({t: tv, Qs: sp.sqrt(q2v)}), 25)
Sp = S_poly(k,3)
Sexpr = sum(sp.Rational(cc.numerator,cc.denominator)*t**ee for ee,cc in enumerate(Sp))
rhs = sp.sqrt(q2)*Sexpr/(2*t*(1-t)*q2**k)
rhs_v = sp.N(rhs.subs(t, tv), 25)
print("lhs(-c/48):", lhs_v)
print("rhs(q S_3/(2t(1-t)q2^4)):", rhs_v)
print("ratio:", sp.N(lhs_v/rhs_v, 10))
# try overall sign and factor variants
print("-rhs:", sp.N(-rhs_v,25))
# maybe claim should have Qs (paper's q) vs q2 powers differently; print raw c
print("raw c:", sp.factor(c))
# direct: lhs expressed as qfree+Q*qlin, get numeric qfree, qlin by solving with two Q values
import mpmath as mp
mp.mp.dps = 30
tvf = mp.mpf(2)/7
qv = mp.sqrt(tvf**2+(1-tvf)**2)
f = lambda qval: complex(sp.N(lhs.subs({t: tv, Qs: sp.Float(qval)}), 25))
v1 = f(qv); v2 = f(-qv)
qfree_num = (v1+v2)/2; qlin_num = (v1-v2)/(2*qv)
print("numeric qfree:", qfree_num, " numeric qlin:", qlin_num)
Sv = sum(sp.Rational(cc.numerator,cc.denominator)*tv**ee for ee,cc in enumerate(Sp))
target = Sv/(2*tv*(1-tv)*q2v**k)
print("my qlin target:", sp.N(target, 12))
