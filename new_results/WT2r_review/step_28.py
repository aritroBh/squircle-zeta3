code = '''
import sys
sys.path.insert(0,'.')
import sympy as sp
from math import comb, factorial
from fractions import Fraction as Fr

# ---- (C1) finite-difference identity, exact check d<=80 ----
def fd_id(d,m):
    return sum((-1)**s*comb(d,s)*s**m*(d-s)**(d-m) for s in range(d+1))
ok1 = all(fd_id(d,m) == factorial(d)*(-1)**m for d in range(1,81) for m in range(0,d+1))
print("(C1) sum_s (-1)^s C(d,s) s^m (d-s)^(d-m) = d! (-1)^m  exact for d=1..80:", ok1)
# and the Stirling lemma used in the proof: sum_s (-1)^{d-s} C(d,s) s^n = d! S(n,d)
def stirling2(n,k):
    if n==k==0: return 1
    if n<k or k<=0: return 0
    return k*stirling2(n-1,k)+stirling2(n-1,k-1)
ok1b = all(sum((-1)**(d-s)*comb(d,s)*s**n for s in range(d+1)) == factorial(d)*stirling2(n,d)
           for d in range(1,25) for n in range(0,30))
print("(C1b) Stirling lemma sum_s (-1)^{d-s} C(d,s) s^n = d! S(n,d), d<=24,n<30:", ok1b)

# ---- (C2) U_d^{(m)}(tau) = i^d (-1)^m C(d,m), exact d<=40 (sympy) ----
I = sp.I
def stirling2(n, k):
    if n < 0 or k < 0: return 0
    if n == k == 0: return 1
    if k == 0 or n == 0: return 0
    prev=[0]*(k+1); prev[0]=1
    for i in range(1,n+1):
        cur=[0]*(k+1)
        for j in range(1,min(i,k)+1): cur[j]=j*prev[j]+prev[j-1]
        prev=cur
    return prev[k] if k<=n else 0
def w(s,m):
    if s==0: return Fr(1) if m==0 else Fr(0)
    if m<s: return Fr(0)
    return Fr(factorial(s)*stirling2(m,s)*2**m, factorial(m))
tau = (1+I)/2
ok2 = True
for d in range(1,41):
    for m in range(0,d+1):
        U = sum(comb(d,s)*sp.Rational(w(s,m))*sp.Rational(w(d-s,d-m))*(1-tau)**(2*s)*tau**(2*(d-s)) for s in range(d+1))
        pred = I**d * (-1)**m * comb(d,m)
        if sp.simplify(U-pred)!=0: ok2=False
print("(C2) U_d^{(m)}((1+i)/2) = i^d (-1)^m C(d,m)  exact for d=1..40:", ok2)

# key algebraic facts behind C2:
print("      tau^2 =", sp.simplify(tau**2), " (1-tau)^2 =", sp.simplify((1-tau)**2), " t(1-t) at v=i:", sp.simplify(tau*(1-tau)))
'''
open('rev_04.py','w').write(code)
import sys, subprocess
work = '<project>/research/08_theory/WT2r_review'
r = subprocess.run([sys.executable, f'{work}/rev_04.py'], capture_output=True, text=True, cwd=work)
print(r.stdout); print(r.stderr[-500:] if r.returncode else "")
