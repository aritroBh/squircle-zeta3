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
def S2(j,m):
    if m>j or m<0: return 0
    if j==0: return 1 if m==0 else 0
    return sum(((-1)**(m-t))*comb(m,t)*t**j for t in range(m+1))//factorial(m)
def Ek_coeffs(k):
    out={}
    def add(j,ia,lb,c): out[(j,ia,lb)]=sp.expand(out.get((j,ia,lb),sp.Integer(0))+c)
    for i in range(k+1):
        j=k-i
        fi=Fr(1,factorial(i)); fj=Fr(1,factorial(j))
        for ta in range(i+1):
            tb=i-ta; coefA=Fr(comb(i,ta))* (-(n-1))**ta * n**tb
            for m in range(j+1):
                bcm = sp.expand(sp.prod([sp.Rational(1,t+1)*(1-2*n-t) for t in range(m)])) if m>0 else sp.Integer(1)
                cm = bcm*factorial(m)*S2(j,m)
                for tal in range(j+1):
                    cl=Fr(comb(j,tal))*Fr(-1)**tal
                    add(m, ta+tal, tb+(j-tal), sp.Rational(fi)*sp.Rational(fj)*coefA*cm*sp.Rational(cl))
    return out
print("Ek_coeffs(1):")
for kk,v in sorted(Ek_coeffs(1).items()): print("  ",kk, sp.expand(v))
# Now check (2n-1) coefficient etc. Expect (1,1,0):(2n-1),(1,0,1):(1-2n),(0,1,0):-(n-1),(0,0,1):n
'''
open(WD+'/dbg1.py','w').write(dbg)
VENV='<project>/.venv/bin/python'
r=subprocess.run([VENV,'dbg1.py'],capture_output=True,text=True,cwd=WD)
print(r.stdout); print(r.stderr[-1500:] if r.returncode else '')