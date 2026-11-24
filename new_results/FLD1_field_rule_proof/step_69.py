import subprocess, time
WD='<project>/research/08_theory/FLD1_field_rule_proof'
VENV='<project>/.venv/bin/python'
# test just mu_formula(6) timing via exec of step_22 top with KMAX=6
code = r'''
import sys, time
sys.argv=['x','6']
src=open("step_22_musymbolic.py").read(); ns={}
exec(src[:src.index("if __name__")], ns)
t=time.time()
F=ns['mu_formula'](6)
print("mu_formula(6) time", round(time.time()-t,1),"s, monomials:", len(F.terms()))
'''
open(WD+'/tmp_k6.py','w').write(code)
t=time.time()
r = subprocess.run(['nice',VENV,'tmp_k6.py'],capture_output=True,text=True,cwd=WD, timeout=1500)
print("wall time", round(time.time()-t,1),"s")
print(r.stdout); print(r.stderr[-1000:] if r.returncode else '')