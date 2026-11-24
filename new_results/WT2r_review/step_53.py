code = '''
import sys
sys.path.insert(0,'.')
sys.path.insert(0, '<project>/paper/v2/code')
import sympy as sp
import mk_integrands
from genintegrand import S_poly
t = mk_integrands.t          # the same symbol object
Qs = mk_integrands.Qs        # their Q (= q)
q2 = t**2+(1-t)**2

def split_Q(expr):
    expr = sp.together(expr)
    N, D = sp.fraction(expr)
    def red(P):
        P = sp.expand(P)
        qf = sp.Integer(0); ql = sp.Integer(0)
        pw = sp.Wild('pw')
        for term in sp.Add.make_args(P):
            mt = term.match(Qs**pw)
            if mt is not None and pw in mt:
                e_ = int(mt[pw]); rest = sp.expand(term/Qs**e_)
            else:
                e_ = 0; rest = term
            # Q^e reduction
            if e_ >= 0:
                f_, l_ = (q2**(e_//2), 0) if e_%2==0 else (0, q2**((e_-1)//2))
            else:
                m = -e_
                f0, l0 = (q2**(m//2), 0) if m%2==0 else (0, q2**((m-1)//2))
                f_, l_ = f0/q2**m, l0/q2**m
            qf += rest*f_; ql += rest*l_
        return sp.expand(qf), sp.expand(ql)
    N0, N1 = red(N); D0, D1 = red(D)
    den2 = sp.expand(D0**2 - D1**2*q2)
    qfree = sp.together(sp.expand(N0*D0 - N1*D1*q2)/den2)
    qlin  = sp.together(sp.expand(N1*D0 - N0*D1)/den2)
    return qfree, qlin

for k in [4,5,6,7]:
    ti = mk_integrands.t_integrand(k)
    bad=[]
    for (i,j),c in ti.items():
        lhs = -c/(2*sp.factorial(k))
        qfree, qlin = split_Q(lhs)
        qf_target = sp.Integer(0)
        if j==0: qf_target += -sp.Rational(1,2)/sp.factorial(k)/t
        if i==0: qf_target += -sp.Rational(1,2)/sp.factorial(k)/(1-t)
        Sp = S_poly(k,i)
        Sexpr = sum(sp.Rational(cc.numerator,cc.denominator)*t**ee for ee,cc in enumerate(Sp))
        ok1 = sp.simplify(qfree - qf_target)==0
        ok2 = sp.simplify(qlin - Sexpr/(2*t*(1-t)*q2**k))==0
        if not (ok1 and ok2): bad.append(((i,j), ok1, ok2))
    print(f"k={k}: genintegrand closed form == paper t_integrand (all {len(ti)} monomials):",
          "ALL MATCH" if not bad else bad[:4], flush=True)
'''
open('rev_13.py','w').write(code)
import sys, subprocess
work = '<project>/research/08_theory/WT2r_review'
r = subprocess.run([sys.executable, f'{work}/rev_13.py'], capture_output=True, text=True, cwd=work, timeout=900)
print(r.stdout); print(r.stderr[-600:] if r.returncode else "")
