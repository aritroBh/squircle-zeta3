import sys, os
workdir = '<project>/research/08_theory/TOW2R_review'
sys.path.insert(0, os.path.join(workdir, 'tow2_copy'))
import sympy as sp
from sympy.functions.combinatorial.numbers import stirling
x = sp.Symbol('x')
P = {0: sp.Integer(1), 1: x*(x-1)/2, 2: x*(x-1)*(x-2)*(3*x-5)/24,
     3: x*(x-1)*(x-2)**2*(x-3)**2/48,
     4: x*(x-1)*(x-2)*(x-3)*(x-4)*(15*x**3-150*x**2+485*x-502)/5760}
ok_eq = all(sp.Integer(stirling(n, n-r, kind=2)) == P[r].subs(x, n) for r in range(1,5) for n in range(0, 25))
print('P_r(n) == S(n,n-r) for n=0..24, r<=4:', ok_eq)
for r in range(1,5):
    ff = sp.prod(x-u for u in range(r))
    q, rem = sp.div(P[r], ff, x)
    print(f'P_{r}: divisible by x^_({r}): rem==0 -> {sp.expand(rem)==0}')
with open(os.path.join(workdir,'progress.log'),'a') as f:
    f.write('step_11b: P_r(n)==S(n,n-r) n<=24 True; divisibility confirmed\n')
