"""Independent P5 derivation and exact direct fifth derivatives at fresh keys."""
import ast, pickle
from math import factorial
from pathlib import Path
import sympy as s
from sympy.functions.combinatorial.numbers import stirling
x,t=s.symbols('x t'); I=s.I
def definitions(path):
    return {n.name:ast.dump(n, include_attributes=False) for n in ast.parse(Path(path).read_text()).body if isinstance(n,ast.FunctionDef)}
a,b=definitions('prove_S34.py'),definitions('prove_S5.py')
for name in ['ff','tf','tf1','normalized_summand','Pi','total']:
    assert a[name]==b[name],name
print('Six finite-sum method functions AST-identical to prove_S34.',flush=True)
# Partitions of the five excess elements in nonsingleton blocks.
P=0
for c2 in range(6):
 for c3 in range(3):
  for c4 in range(2):
   for c5 in range(2):
    for c6 in range(2):
     cs=[c2,c3,c4,c5,c6]
     if sum((j-1)*c for j,c in zip(range(2,7),cs))!=5: continue
     used=sum(j*c for j,c in zip(range(2,7),cs))
     P+=s.prod(x-u for u in range(used))/s.prod(factorial(j)**c*factorial(c) for j,c in zip(range(2,7),cs))
expected=x*(x-1)*(x-2)*(x-3)*(x-4)**2*(x-5)**2*(3*x*x-23*x+38)/11520
assert s.expand(P-expected)==0
for n in range(6,36): assert P.subs(x,n)==stirling(n,n-5,kind=2)
print('Independent block-count P5 polynomial matches; Stirling integer checks n=6..35 pass.',flush=True)
raw=pickle.load(open('../TOW5_tower_depth10/N5_proved.pkl','rb'))
d,m=s.symbols('d m'); N=s.sympify(raw['N5'])
T=s.Poly(t,t,domain=s.QQ_I); U=1-T; W=T*T+U*U
for dd,pp in [(11,11),(11,7),(13,8),(14,7),(17,12)]:
 qq=dd-pp; out=s.Poly(0,t,domain=s.QQ_I)
 for j in range(1,dd+1):
  C=s.binomial(s.Rational(1,2),j)*factorial(j)*2**dd
  for r in range(max(0,j-qq),min(j,pp)+1):
   c=C*stirling(pp,r,kind=2)*stirling(qq,j-r,kind=2)
   if c: out+=c*T**(2*(j-r))*U**(2*r)*W**(dd-j)
 ts=(1+I)/2
 direct=s.cancel(out.diff((t,5)).eval(ts)/out.eval(ts))
 pred=s.cancel(N.subs({d:dd,m:pp})/(3*s.prod(2*dd-2*j-1 for j in range(1,6))))
 assert s.expand(direct-pred)==0,(dd,pp)
 print('Fresh direct exact differentiation PASS',dd,pp,flush=True)
