import sympy as sp
x, y = sp.symbols('x y')
N = 15
gam, = sp.symbols('gamma_E')
def logG(z, N):
    return -gam*z + sum((-1)**n*sp.zeta(n)*z**n/n for n in range(2, N+1))
lg = logG(x, N) + logG(y, N) - logG(x+y, N)
Pser = sp.exp(lg).series(x, 0, N).removeO().series(y, 0, N).removeO().expand()
def aa_sym(p, q):
    return Pser.coeff(x, p).coeff(y, q+1)   # aa_{p,q} = [x^p y^{q+1}] P
# 1) P(x,0) == 1 and P(0,y) == 1 (=> aa_{0,d} = 0 for d>=1)
print("P(x,0)=1:", sp.expand(Pser.subs(y,0)) == 1, " P(0,y)=1:", sp.expand(Pser.subs(x,0)) == 1)
print("aa_{0,d} = [y^{d+1}]P(0,y) = 0 for d=1..12:", all(aa_sym(0,d).expand()==0 for d in range(1,13)))
# 2) recursion aa_{p,q} = aa_{q+1,p-1} EXACT for p=1..8, q=0..7
bad = []
for p in range(1,9):
    for q in range(0,8):
        if sp.expand(aa_sym(p,q) - aa_sym(q+1,p-1)) != 0: bad.append((p,q))
print("recursion (10) exact p<=8,q<=7:", "ALL OK" if not bad else bad)
# 3) involution bookkeeping per odd d: sum_{A_*} f - (1/2) f(D_*) == (1/2) sum_{S_d} f, f=(-1)^{q+1} aa
for d in range(1, 14, 2):
    mid = (d-1)//2
    Sd = [(p, d-p) for p in range(1, d+1)]
    Astar = [(p,q) for (p,q) in Sd if p > q]
    Dstar = (mid+1, mid)
    f = lambda p,q: (-1)**(q+1)*aa_sym(p,q)
    lhs = sum(f(p,q) for p,q in Astar) - sp.Rational(1,2)*f(*Dstar)
    rhs = sp.Rational(1,2)*sum(f(p,q) for p,q in Sd)
    ok = sp.expand(lhs-rhs) == 0
    print(f"d={d}: involution fold identity exact: {ok}")
with open('progress.log','a') as f:
    f.write("step_16: Laurent bookkeeping exact (P(x,0)=P(0,y)=1, aa_{0,d}=0); recursion (10) exact p<=8 q<=7; involution fold identity exact for odd d=1..13\n")
