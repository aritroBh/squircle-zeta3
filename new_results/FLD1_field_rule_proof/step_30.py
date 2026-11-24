import subprocess
WD='<project>/research/08_theory/FLD1_field_rule_proof'
# import the module functions from step_22 by exec, then print mu_formula(1)
dbg = r'''
import sys
sys.path.insert(0,'.')
import importlib.util
spec = importlib.util.spec_from_file_location("ms", "step_22_musymbolic.py")
# We cannot run __main__; instead re-create minimal pieces by exec of the top part
src = open("step_22_musymbolic.py").read()
top = src[:src.index("if __name__")]
ns = {}
exec(top, ns)
F = ns['mu_formula'](1)
print("mu_{1,n} terms:")
for mon, cf in F.terms():
    print("   mon", mon, " coef", sp.factor(ns['K'].to_sympy(cf) if hasattr(ns['K'],'to_sympy') else cf.as_expr()))
import sympy as sp
'''
open(WD+'/dbg5.py','w').write(dbg)
VENV='<project>/.venv/bin/python'
r=subprocess.run([VENV,'dbg5.py'],capture_output=True,text=True,cwd=WD)
print(r.stdout); print(r.stderr[-1500:] if r.returncode else '')