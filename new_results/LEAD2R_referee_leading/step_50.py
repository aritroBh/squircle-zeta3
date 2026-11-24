import sympy as sp
from math import factorial
N = 14
gam = sp.symbols('gE')
# lg as dict of (i,j) -> coeff, total degree <= N
lg = {}
for n in range(1, N+1):
    c = -gam if n == 1 else (-1)**n*sp.zeta(n)/n
    # logG(x): (n,0); logG(y): (0,n); -logG(x+y): -(n-k,k) binom(n,k)
    lg[(n,0)] = lg.get((n,0),0) + c
    lg[(0,n)] = lg.get((0,n),0) + c
    for k in range(n+1):
        lg[(n-k,k)] = lg.get((n-k,k),0) - c*sp.binomial(n,k)
lg = {k: sp.expand(v) for k, v in lg.items() if k[0]+k[1] <= N}
def polymul(A, B, N):
    C = {}
    for (i,j), a in A.items():
        for (k,l), b in B.items():
            if i+k+j+l <= N:
                C[(i+k,j+l)] = C.get((i+k,j+l),0) + a*b
    return C
# P = exp(lg) truncated
P = {(0,0): sp.Integer(1)}
term = {(0,0): sp.Integer(1)}
for m in range(1, N+1):
    term = polymul(term, lg, N)
    for k, v in term.items():
        P[k] = P.get(k,0) + v/factorial(m)
P = {k: sp.expand(v) for k, v in P.items()}
def aa_sym(p,q): return P.get((p,q+1), sp.Integer(0))
print("P(x,0)=1:", all(sp.expand(P.get((i,0),0))==0 for i in range(1,N+1)))
print("P(0,y)=1:", all(sp.expand(P.get((0,j),0))==0 for j in range(1,N+1)))
print("aa_{0,d}=0, d=1..12:", all(sp.expand(aa_sym(0,d))==0 for d in range(1,13)))
bad = [(p,q) for p in range(1,9) for q in range(0,8) if sp.expand(aa_sym(p,q)-aa_sym(q+1,p-1))!=0]
print("recursion (10) exact p<=8,q<=7:", "ALL OK" if not bad else bad)
res = {}
for d in range(1, 14, 2):
    mid=(d-1)//2
    Sd=[(p,d-p) for p in range(1,d+1)]
    Astar=[(p,q) for p,q in Sd if p>q]
    Dstar=(mid+1,mid)
    f=lambda p,q: (-1)**(q+1)*aa_sym(p,q)
    lhs = sum(f(p,q) for p,q in Astar) - sp.Rational(1,2)*f(*Dstar)
    rhs = sp.Rational(1,2)*sum(f(p,q) for p,q in Sd)
    res[d] = (sp.expand(lhs-rhs)==0, sp.expand(rhs))
    print(f"d={d}: fold identity exact: {res[d][0]},  (1/2)sum_Sd f = {res[d][1]}")
with open('progress.log','a') as f:
    f.write("step_17: exact symbolic (dict-poly) aa_{p,q}; P(x,0)=P(0,y)=1; aa_{0,d}=0; recursion (10) exact p<=8 q<=7; fold identity exact d=1,3,..,13\n")
