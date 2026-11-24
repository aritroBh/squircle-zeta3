
import sys
sys.path.insert(0,'.')
sys.path.insert(0, '<project>/paper/v2/code')
import sympy as sp
import mk_integrands
from genintegrand import S_poly
from fractions import Fraction as Fr
t, Q = sp.symbols('t Q')
q2 = t**2+(1-t)**2
def split_Q(e):
    e = sp.expand(e)
    qf = sp.Integer(0); ql = sp.Integer(0)
    for term in sp.Add.make_args(e):
        if not term.has(Q): qf += term; continue
        p = sp.Poly(term, Q)
        for mono, cf in p.terms():
            deg = mono[0]
            if deg % 2 == 0: qf += cf*q2**(deg//2)
            else: ql += cf*q2**((deg-1)//2)
    return qf, ql
k=4; key=(3,1)
ti = mk_integrands.t_integrand(k)
c = ti[key]
print("raw paper coeff:", sp.factor(c))
lhs = -c/(2*sp.factorial(k))
qfree, qlin = split_Q(lhs)
tv = Fr(2,7)
def ev(e):
    e = sp.expand(e)
    val = e.subs(t, sp.Rational(tv.numerator, tv.denominator))
    return sp.nsimplify(val, rational=True)
q2v = tv**2+(1-tv)**2
Sp = S_poly(k,key[0])
Sv = sum(Fr(cc.numerator,cc.denominator)*tv**ee for ee,cc in enumerate(Sp))
print("paper qfree at t=2/7:", ev(qfree))
print("my A/B-delta target at t=2/7:", 0)
print("paper qlin at t=2/7:", ev(qlin))
mine = Sv/(2*tv*(1-tv)*q2v**k)
print("my qlin target at t=2/7:", sp.Rational(mine.numerator, mine.denominator))
print("ratio paper/mine:", sp.simplify(ev(qlin)/sp.Rational(mine.numerator, mine.denominator)))
