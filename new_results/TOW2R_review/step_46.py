import sys, os
workdir = '<project>/research/08_theory/TOW2R_review'
sys.path.insert(0, os.path.join(workdir, 'tow2_copy'))
import sympy as sp
x = sp.Symbol('x')
P = {0: sp.Integer(1), 1: x*(x-1)/2, 2: x*(x-1)*(x-2)*(3*x-5)/24,
     3: x*(x-1)*(x-2)**2*(x-3)**2/48,
     4: x*(x-1)*(x-2)*(x-3)*(x-4)*(15*x**3-150*x**2+485*x-502)/5760}
# 1. recurrence P_r(x+1) = (x+1-r) P_{r-1}(x) + P_r(x), symbolic
ok_rec = all(sp.expand(P[r].subs(x,x+1) - ((x+1-r)*P[r-1] + P[r])) == 0 for r in range(1,5))
print('recurrence r<=4 symbolic:', ok_rec)
# 2. anchors P_r(0..r) hmm: P_r(m) = S(m,m-r): anchors P_r(0)=...=P_r(r-1)=0 (S(m,neg)=0... S(m,m-r) for m<r: m-r<0 -> 0) and P_r(r)=1 (S(r,0)=... S(r,0)=0 for r>=1! hmm)
# Careful: S(m, m-r) at m=r: S(r,0) = 0 (r>=1). At m=r+1: S(r+1,1) = 1.
ok_anch = all(P[r].subs(x, mm) == 0 for r in range(1,5) for mm in range(0, r+1))
print('anchors P_r(m)=0 for 0<=m<=r:', ok_anch)
# 3. equality with Stirling S(n, n-r) at integers n = 0..24
from sympy import stirling
ok_eq = all(sp.Integer(stirling(n, n-r, kind=2)) == P[r].subs(x, n) for r in range(1,5) for n in range(0, 25))
print('P_r(n) == S(n,n-r) for n=0..24, r<=4:', ok_eq)
# 4. divisibility by falling factorial x^underline(r) (needed for removable singularity)
for r in range(1,5):
    ff = sp.prod(x-u for u in range(r))
    q, rem = sp.div(P[r], ff, x)
    print(f'P_{r} / x^underline({r}): remainder {sp.expand(rem)}, quotient poly: {q}')
with open(os.path.join(workdir,'progress.log'),'a') as f:
    f.write('step_11: Stirling tails P_r (r<=4): recurrence+anchors+equality at n<=24 verified; divisible by x^_(r)\n')
