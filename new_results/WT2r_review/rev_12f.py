
import sys
sys.path.insert(0,'.')
sys.path.insert(0, '<project>/paper/v2/code')
import sympy as sp
import mk_integrands
from genintegrand import S_poly
t, Q = sp.symbols('t Q')
q2 = t**2+(1-t)**2

def red_pow(e):
    """Q^e (e integer, possibly negative) -> (qfree_part, qlin_part) polys in t"""
    if e >= 0:
        if e % 2 == 0: return q2**(e//2), 0
        return 0, q2**((e-1)//2)
    m = -e
    # 1/Q^m = Q^{m}/q2^{m}  (since Q^{2m}=q2^m => Q^{-m} = Q^m/q2^m)
    return red_pow(m)[0]/q2**m, red_pow(m)[1]/q2**m

def split_Q(expr):
    expr = sp.together(expr)
    N, D = sp.fraction(expr)
    def red(P):
        P = sp.expand(P)
        qf = sp.Integer(0); ql = sp.Integer(0)
        for term in sp.Add.make_args(P):
            pw = sp.Wild('pw')
            mt = term.match(Q**pw)
            if mt is not None and pw in mt:
                e_ = int(mt[pw]); rest = sp.simplify(term/Q**e_)
            else:
                e_ = 0; rest = term
            f_, l_ = red_pow(e_)
            qf += rest*f_; ql += rest*l_
        return sp.expand(qf), sp.expand(ql)
    N0, N1 = red(N); D0, D1 = red(D)
    den2 = sp.expand(D0**2 - D1**2*q2)
    qfree = sp.expand(N0*D0 - N1*D1*q2)/den2
    qlin  = sp.expand(N1*D0 - N0*D1)/den2
    return sp.together(qfree), sp.together(qlin)

# sanity on the k=4 (3,1) raw coefficient from before
ti4 = mk_integrands.t_integrand(4)
c = ti4[(3,1)]
qf, ql = split_Q(-c/(2*sp.factorial(4)))
Sp = S_poly(4,3)
Sexpr = sum(sp.Rational(cc.numerator,cc.denominator)*t**ee for ee,cc in enumerate(Sp))
print("sanity k=4 (3,1): qfree==0:", sp.simplify(qf)==0,
      " qlin==S_3/(2t(1-t)q2^4):", sp.simplify(ql - Sexpr/(2*t*(1-t)*q2**4))==0)
