code = '''
import sys, time
sys.path.insert(0,'.')
sys.path.insert(0, '<project>/paper/v2/code')
import sympy as sp
import mk_integrands
from genintegrand import S_poly
t, Q = sp.symbols('t Q', positive=True)
k = 6
t0=time.time()
ti = mk_integrands.t_integrand(k)
print(f"t_integrand({k}) computed in {time.time()-t0:.1f}s, {len(ti)} terms", flush=True)
q2 = t**2+(1-t)**2
bad=[]
for (i,j),c in ti.items():
    # M_k coefficient = -c/(2 k!) ; claim sqrt part = q S_i/(2t(1-t)q2^k) with Q->q, plus A/B deltas
    claim = Q*S_poly(k,i)/(2*t*(1-t)*q2**k)
    if j==0: claim += -1/(2*sp.factorial(k)*t)
    if i==0: claim += -1/(2*sp.factorial(k)*(1-t))
    lhs = sp.simplify(-c/(2*sp.factorial(k)))
    rhs = sp.simplify(claim.subs(Q, sp.sqrt(q2)))
    if sp.simplify(lhs-rhs)!=0: bad.append(((i,j), lhs, rhs))
print("genintegrand closed form == paper t_integrand(k=6) all monomials:", "ALL MATCH" if not bad else bad[:3])
'''
open('rev_12.py','w').write(code)
import sys, subprocess
work = '<project>/research/08_theory/WT2r_review'
r = subprocess.run([sys.executable, f'{work}/rev_12.py'], capture_output=True, text=True, cwd=work, timeout=900)
print(r.stdout); print(r.stderr[-600:] if r.returncode else "")
