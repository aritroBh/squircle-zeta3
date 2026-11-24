import subprocess
WD='<project>/research/08_theory/FLD1_field_rule_proof'
dbg = r'''
import sys
sys.path.insert(0,'.')
from sympy.polys.fields import field
from sympy.polys.rings import ring
from math import comb, factorial
import zsum
from zsum import K, nn
KMAX=3
Hv=[f'h{r}' for r in range(1,KMAX+1)]; Gv=[f'g{r}' for r in range(1,KMAX+1)]; Zv=[f'z{r}' for r in range(2,KMAX+1)]
R=ring(', '.join(Hv+Gv+Zv), K)[0]
gens=R.gens
h={r:gens[r-1] for r in range(1,KMAX+1)}
g={r:gens[KMAX+r-1] for r in range(1,KMAX+1)}
z={r:gens[2*KMAX+r-2] for r in range(2,KMAX+1)}
def shift_h(r,j): return h[r]+sum((K(1)/(nn+t)**r for t in range(j)),K(0))
def shift_g(r,j): return g[r]+sum((K(1)/(2*nn+t)**r for t in range(j)),K(0))
def psi(m,which,j):
    sh=shift_h if which=='h' else shift_g
    if m==0: return sh(1,j)
    c=(-1)**(m+1)*factorial(m)
    harm=sh(m+1,j) if (m+1)<=KMAX else None
    out=K(0)
    if (m+1) in z: out+=c*z[m+1]
    if harm is not None: out+=-c*harm
    return out
print("psi(0,'h',1)=", psi(0,'h',1))
print("psi(0,'g',1)=", psi(0,'g',1))
print("cumulant(1,0,1)=", -psi(0,'g',1)+psi(0,'h',1))
# now the moments version
def cumulant(r,m,j):
    out=-psi(r+m-1,'g',j)
    if m==0: out=out+psi(r+m-1,'h',j)
    if r==0: out=out+psi(r+m-1,'h',0)
    return out
print("cumulant(1,0,1) via func =", cumulant(1,0,1))
'''
open(WD+'/dbg4.py','w').write(dbg)
VENV='<project>/.venv/bin/python'
r=subprocess.run([VENV,'dbg4.py'],capture_output=True,text=True,cwd=WD)
print(r.stdout); print(r.stderr[-1500:] if r.returncode else '')