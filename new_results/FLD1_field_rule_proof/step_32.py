import subprocess
WD='<project>/research/08_theory/FLD1_field_rule_proof'
dbg = r'''
import sys
sys.path.insert(0,'.')
src = open("step_22_musymbolic.py").read()
top = src[:src.index("if __name__")]
ns = {}
exec(top, ns)
sp=ns['sp']; K=ns['K']; nn=ns['nn']; R=ns['R']
Ek=ns['Ek_coeffs']; moments=ns['moments']; toK=ns['toK']
k=1
Dc={}
tot=R(0)
for (j,ia,lb),co in Ek(k).items():
    if j not in Dc: Dc[j]=moments(j)
    pr=K(1)
    for t in range(j): pr*=(nn+t)/(2*nn+t)
    D=Dc[j][(ia,lb)]
    term = toK(co)*pr*D
    print(f"(j={j},ia={ia},lb={lb}) co={sp.expand(co)} pr={pr.as_expr()} D={D}")
    print("    toK(co)=", toK(co).as_expr(), " term=", term)
    tot+=term
print("TOT (before prefactor) =", tot)
print("prefactor -1/(n(2n-1)) =", (-K(1)/(nn*(2*nn-1))).as_expr())
print("mu =", tot*(-K(1)/(nn*(2*nn-1))))
'''
open(WD+'/dbg7.py','w').write(dbg)
VENV='<project>/.venv/bin/python'
r=subprocess.run([VENV,'dbg7.py'],capture_output=True,text=True,cwd=WD)
print(r.stdout); print(r.stderr[-1500:] if r.returncode else '')