import subprocess
WD='<project>/research/08_theory/FLD1_field_rule_proof'
dbg = r'''
import sys
sys.path.insert(0,'.')
import sympy as sp
from sympy import QQ
from sympy.polys.fields import field
from sympy.polys.rings import ring
from fractions import Fraction as Fr
from math import comb, factorial
import zsum
from zsum import K, nn, ONE, ZERO
KMAX=3
n=sp.Symbol('n')
Hv=[f'h{r}' for r in range(1,KMAX+1)]; Gv=[f'g{r}' for r in range(1,KMAX+1)]; Zv=[f'z{r}' for r in range(2,KMAX+1)]
R=ring(', '.join(Hv+Gv+Zv), K)[0]
gens=R.gens
h={r:gens[r-1] for r in range(1,KMAX+1)}
g={r:gens[KMAX+r-1] for r in range(1,KMAX+1)}
z={r:gens[2*KMAX+r-2] for r in range(2,KMAX+1)}
def shift_h(r,j): return h[r]+sum((K(1)/(nn+t)**r for t in range(j)),K(0))
def shift_g(r,j): return g[r]+sum((K(1)/(2*nn+t)**r for t in range(j)),K(0))
def psi(m,which,j):
    H=h if which=='h' else g
    sh=shift_h if which=='h' else shift_g
    if m==0: return sh(1,j)
    c=(-1)**(m+1)*factorial(m)
    harm=sh(m+1,j) if (m+1)<=KMAX else None
    out=K(0)
    if (m+1) in z: out+=c*z[m+1]
    if harm is not None: out+=-c*harm
    return out
def cumulant(r,m,j):
    out=-psi(r+m-1,'g',j)
    if m==0: out=out+psi(r+m-1,'h',j)
    if r==0: out=out+psi(r+m-1,'h',0)
    return out
def moments(j):
    D={(0,0):R(1)}
    for tot in range(1,KMAX+1):
        for I in range(tot,-1,-1):
            L=tot-I; v=R(0)
            if I>0:
                for r in range(1,I+1):
                    for m in range(0,L+1):
                        v+=comb(I-1,r-1)*comb(L,m)*cumulant(r,m,j)*D[(I-r,L-m)]
            else:
                for m in range(1,L+1):
                    v+=comb(L-1,m-1)*cumulant(0,m,j)*D[(0,L-m)]
            D[(I,L)]=v
    return D
D0=moments(0)
print("D_{1,0}(n,n) =", D0[(1,0)])
print("D_{0,1}(n,n) =", D0[(0,1)])
D1=moments(1)
print("D_{1,0}(n+1,n) =", D1[(1,0)])
print("D_{0,1}(n+1,n) =", D1[(0,1)])
'''
open(WD+'/dbg2.py','w').write(dbg)
VENV='<project>/.venv/bin/python'
r=subprocess.run([VENV,'dbg2.py'],capture_output=True,text=True,cwd=WD)
print(r.stdout); print(r.stderr[-1500:] if r.returncode else '')