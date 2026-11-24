code = '''
import sys
sys.path.insert(0,'.')
sys.path.insert(0, '<project>/paper/v2/code')
import sympy as sp
import mk_integrands
from genintegrand import S_poly
t, Q = sp.symbols('t Q')
q2 = t**2+(1-t)**2

def split_Q(e):
    """expr (rational, Q algebraic with Q^2=q2) -> (qfree, qlin) exact"""
    e = sp.expand(e)
    terms = sp.Add.make_args(e)
    qf = sp.Integer(0); ql = sp.Integer(0)
    for term in terms:
        # extract power of Q
        p = sp.Poly(term, Q) if term.has(Q) else None
        if p is None:
            qf += term; continue
        degs = [m[0] for m in p.terms()]
        for deg in degs:
            cf = p.coeff_monomial(deg)
            if deg % 2 == 0:
                qf += cf * q2**(deg//2)
            else:
                ql += cf * q2**((deg-1)//2)
    return sp.expand(qf), sp.expand(ql)

for k in [4,5,6]:
    ti = mk_integrands.t_integrand(k)
    bad=[]
    for (i,j),c in ti.items():
        lhs = sp.expand(-c/(2*sp.factorial(k)))
        qfree, qlin = split_Q(lhs)
        qf_target = sp.Integer(0)
        if j==0: qf_target += -sp.Rational(1,2)/sp.factorial(k)/t
        if i==0: qf_target += -sp.Rational(1,2)/sp.factorial(k)/(1-t)
        Sp = S_poly(k,i)
        Sexpr = sum(sp.Rational(cc.numerator,cc.denominator)*t**ee for ee,cc in enumerate(Sp))
        ok1 = sp.simplify(sp.together(qfree - qf_target))==0
        ok2 = sp.simplify(sp.together(qlin - Sexpr/(2*t*(1-t)*q2**k)))==0
        if not (ok1 and ok2): bad.append(((i,j), ok1, ok2))
    print(f"k={k}: genintegrand == paper t_integrand:", "ALL MATCH" if not bad else bad[:4], flush=True)
'''
open('rev_12d.py','w').write(code)
import sys, subprocess
work = '<project>/research/08_theory/WT2r_review'
r = subprocess.run([sys.executable, f'{work}/rev_12d.py'], capture_output=True, text=True, cwd=work, timeout=600)
print(r.stdout); print(r.stderr[-600:] if r.returncode else "")
