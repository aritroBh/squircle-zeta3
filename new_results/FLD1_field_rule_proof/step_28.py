import subprocess
WD='<project>/research/08_theory/FLD1_field_rule_proof'
dbg = r'''
import sys
sys.path.insert(0,'.')
from sympy.polys.fields import field
from sympy.polys.rings import ring
import zsum
from zsum import K, nn
KMAX=3
Hv=[f'h{r}' for r in range(1,KMAX+1)]; Gv=[f'g{r}' for r in range(1,KMAX+1)]; Zv=[f'z{r}' for r in range(2,KMAX+1)]
R=ring(', '.join(Hv+Gv+Zv), K)[0]
gens=R.gens
h={r:gens[r-1] for r in range(1,KMAX+1)}
g={r:gens[KMAX+r-1] for r in range(1,KMAX+1)}
print("shift_h(1,1) =", h[1] + sum((K(1)/(nn+t)**1 for t in range(1)), K(0)))
print("shift_g(1,1) =", g[1] + sum((K(1)/(2*nn+t)**1 for t in range(1)), K(0)))
print("h[1]=", h[1], " g[1]=", g[1])
print("list gens:", list(gens))
'''
open(WD+'/dbg3.py','w').write(dbg)
VENV='<project>/.venv/bin/python'
r=subprocess.run([VENV,'dbg3.py'],capture_output=True,text=True,cwd=WD)
print(r.stdout); print(r.stderr[-1500:] if r.returncode else '')