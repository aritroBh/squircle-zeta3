import subprocess
WD='<project>/research/08_theory/FLD1_field_rule_proof'
dbg = r'''
import sys
sys.path.insert(0,'.')
import sympy as sp
src = open("step_22_musymbolic.py").read()
top = src[:src.index("if __name__")]
ns = {}
exec(top, ns)
sp = ns['sp']
for kk in (1,2):
    F = ns['mu_formula'](kk)
    print(f"mu_{{{kk},n}} terms:")
    for mon, cf in F.terms():
        print("   mon", mon, " coef", sp.factor(cf.as_expr()))
    # evaluate at n=1 with h,g harmonic
    from fractions import Fraction as Fr
    def Hf(r,m): return sum((Fr(1,i**r) for i in range(1,m+1)),Fr(0))
    KMAX=ns['KMAX']
    def evalK(c,nv): return Fr(int(c.numer(nv)),1)/Fr(int(c.denom(nv)),1)
    for nv in (1,2):
        acc=Fr(0)
        for mon,cf in F.terms():
            cv=evalK(cf,nv); hv=Fr(1)
            for r in range(1,KMAX+1):
                hv*=Hf(r,nv-1)**mon[r-1]; hv*=Hf(r,2*nv-1)**mon[KMAX+r-1]
            acc+=cv*hv
        print(f"   mu_{{{kk},{nv}}} (rational+zeta-free part if no z) =", acc)
'''
open(WD+'/dbg6.py','w').write(dbg)
VENV='<project>/.venv/bin/python'
r=subprocess.run([VENV,'dbg6.py'],capture_output=True,text=True,cwd=WD)
print(r.stdout); print(r.stderr[-1500:] if r.returncode else '')