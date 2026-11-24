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
tv = sp.Rational(2,7)
for key in [(3,1),(4,0),(2,2)]:
    c = ti4[key]
    lhs = -c/(2*sp.factorial(4))
    lhs_v = complex(lhs.subs({t: tv, Q: sp.sqrt(q2)}).evalf(30))
    for mname, m in [("S_i (i=key[0])", key[0]), ("S_j", key[1])]:
        Sp = S_poly(4,m)
        Sexpr = sum(sp.Rational(cc.numerator,cc.denominator)*t**ee for ee,cc in enumerate(Sp))
        rhs_v = complex((sp.sqrt(q2)*Sexpr/(2*t*(1-t)*q2**4)).subs(t, tv).evalf(30))
        print(f"key={key} {mname}: lhs={lhs_v:.12g} rhs={rhs_v:.12g} ratio={lhs_v/rhs_v:.6g}")
    # A/B delta for (4,0),(0,4):
    print()
'''
open('rev_12h.py','w').write(code)
import sys, subprocess
work = '<project>/research/08_theory/WT2r_review'
r = subprocess.run([sys.executable, f'{work}/rev_12h.py'], capture_output=True, text=True, cwd=work, timeout=300)
print(r.stdout); print(r.stderr[-600:] if r.returncode else "")
