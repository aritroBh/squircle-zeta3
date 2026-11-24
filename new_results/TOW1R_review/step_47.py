import subprocess, os
work = "<project>/research/08_theory/TOW1R_review/"
step = r'''
import sys, os
sys.dont_write_bytecode=True
WTS="<project>/research/08_theory/WTS_weight_k_structure/theory/"
sys.path.insert(0,WTS)
from fractions import Fraction as Fr
from genintegrand import S_poly
import sympy as sp
i=sp.I; v=sp.symbols('v'); sq2=sp.sqrt(2); a=1+sq2; b=sq2-1
def A_sym(d,m):
    P=S_poly(d,m); D=len(P)-1
    tv=(v**2+2*v-1)/(4*v)
    St=sum(Fr(c)*tv**s for s,c in enumerate(P))*(4*v)**D
    return -(sq2/2)*8**d/4**D * v**(2*d-1-D) * St * (v**2+1)**(2-2*d)/((v**2-a**2)*(v**2-b**2))
bad1=bad2=0
for d in range(2,10):
    for m in range(1,d+1):
        A=A_sym(d,m)
        A1=sp.diff(A,v); A2=sp.diff(A,v,2)
        # A'(i)/A(i) via simplify of rational ratio
        r1 = sp.simplify(sp.expand(A1)/sp.expand(A))
        val1 = sp.simplify(r1.subs(v,i))
        if sp.simplify(val1 - (-i*(d-2)))!=0: bad1+=1; print("(logA)' mismatch",d,m,val1)
        r2 = sp.simplify(sp.expand(A2)/sp.expand(A))
        val2 = sp.simplify(r2.subs(v,i))
        target = sp.simplify((3*d-3)/2 - i/2*(2*(d-2*m)-2*d*i-4*(d-1)*(d-2*m)/(2*d-3)) - (d-2)**2)
        if sp.simplify(val2-target)!=0: bad2+=1; print("A''/A mismatch",d,m,sp.simplify(val2),target)
print("(log A)'(i) = -i(d-2), d=2..9:", "0 bad" if bad1==0 else f"{bad1} bad")
print("A''/A(i) = 2 rho2/rho0, d=2..9:", "0 bad" if bad2==0 else f"{bad2} bad")
'''
with open(work+"step_19.py","w") as f: f.write(step)
env=dict(os.environ); env["PYTHONDONTWRITEBYTECODE"]="1"
r=subprocess.run(["nice","-n","19","<project>/.venv/bin/python","step_19.py"],
                 cwd=work,capture_output=True,text=True,env=env,timeout=900)
print(r.stdout); print("STDERR:",r.stderr[:2500])
with open(work+"progress.log","a") as fh: fh.write("step_19: rho1/rho2 derivation via direct A'/A, A''/A: verified d=2..9\n")
