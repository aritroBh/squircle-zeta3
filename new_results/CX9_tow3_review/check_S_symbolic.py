"""Independent h=t-t* coefficient derivation of R1..R4."""
from pathlib import Path
from math import comb,factorial
import sys
import sympy as sp
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'TOW3_tower_depth6'))
import tow3lib as T
d,m=T.d,T.m; i=sp.I; x=sp.Symbol('x')
ts=(1+i)/2; zs=1-ts
P={0:sp.Integer(1),1:x*(x-1)/2,2:x*(x-1)*(x-2)*(3*x-5)/24,
3:x*(x-1)*(x-2)**2*(x-3)**2/48,
4:x*(x-1)*(x-2)*(x-3)*(x-4)*(15*x**3-150*x**2+485*x-502)/5760}
for r in range(1,5):
    assert sp.expand(P[r].subs(x,x+1)-(x+1-r)*P[r-1]-P[r])==0
    assert all(P[r].subs(x,a)==0 for a in range(r+1))
print('PASS Stirling-tail recurrence and anchors through r=4',flush=True)
# q2(ts+h)=2*i*h+2*h^2. Cancel binomial/falling-factorial
# ratios before evaluating endpoints m=0,d; no 0/0 boundary substitutions.
for l in range(1,5):
    result=0
    for n in range(l+1):
        for r in range(n+1):
            base=(-2)**n/T.Pi(n)*P[r].subs(x,m)*P[n-r].subs(x,d-m)*ts**(2*(r-n))*zs**(-2*r)
            aa,bb=2*(d-m-n+r),2*(m-r)
            for j in range(min(n,l-n)+1):
                rem=l-n-j
                bc=sum(sp.ff(aa,c)/factorial(c)/ts**c*sp.ff(bb,rem-c)/factorial(rem-c)*(-1/zs)**(rem-c) for c in range(rem+1))
                result+=base*2**n*comb(n,j)*i**(n-j)*bc
    assert sp.cancel(sp.expand(factorial(l)*result)-T.T.R[l])==0,l
    print(f'PASS independent symbolic S^{l}/S coefficient derivation in d,m',flush=True)
