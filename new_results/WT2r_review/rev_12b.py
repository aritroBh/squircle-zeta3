
import sys, time
sys.path.insert(0,'.')
sys.path.insert(0, '<project>/paper/v2/code')
import sympy as sp
import mk_integrands
from genintegrand import S_poly
t, Q = sp.symbols('t Q', positive=True)
q2 = t**2+(1-t)**2

def reduce_Q(expr):
    """collect powers of Q, replace Q^2 -> q2, return (qfree, qlin) with expr = qfree + Q*qlin"""
    e = sp.expand(expr)
    e = e.replace(Q**2, q2)  # may need repeated
    for _ in range(12):
        e = sp.expand(e.replace(Q**2, q2))
    qfree = sp.expand(e.subs(Q, 0))
    qlin = sp.expand((e - qfree)/Q) if e != qfree else sp.Integer(0)
    qlin = sp.expand(qlin.subs(Q,0))
    return qfree, qlin

for k in [5,6]:
    ti = mk_integrands.t_integrand(k)
    bad=[]
    for (i,j),c in ti.items():
        m = i
        Sp = S_poly(k,m)
        Sexpr = sum(sp.Rational(cc.numerator,cc.denominator)*t**e for e,cc in enumerate(Sp))
        # claim: M_k coeff = -c/(2k!) = Q*S_m/(2t(1-t)q2^k) + deltas
        qf_target = 0
        if j==0: qf_target += -sp.Rational(1,2)/sp.factorial(k)/t
        if i==0: qf_target += -sp.Rational(1,2)/sp.factorial(k)/(1-t)
        lhs = -c/(2*sp.factorial(k))
        qfree, qlin = reduce_Q(lhs)
        ok1 = sp.simplify(qfree - qf_target)==0
        ok2 = sp.simplify(qlin - Sexpr/(2*t*(1-t)*q2**k))==0
        if not (ok1 and ok2): bad.append(((i,j), ok1, ok2))
    print(f"k={k}: genintegrand == paper t_integrand (Q^2=q2 reduced):", "ALL MATCH" if not bad else bad[:4], flush=True)
