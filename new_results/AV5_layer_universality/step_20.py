import sympy as sp, mpmath, pickle, json
from mpmath import mp, mpf, sqrt, zeta, log, exp, quad, tanh, nstr, factorial, polylog
mp.dps = 90
L2 = log(2)
def eta(n): return (1-mpf(2)**(1-n))*zeta(n)
g1 = lambda s: log(1+exp(-2*s))
g1p = lambda s: -2*exp(-2*s)/(1+exp(-2*s))

def Iab(a,b):
    return quad(lambda s: s**a*g1(s)**b, [0,mpf('0.5'),1,2,4,8,16,40])
def K2int(m):
    return quad(lambda s: s**m*g1(s)*g1p(s)**2, [0,mpf('0.5'),1,2,4,8,16,40])
def K3int(m):
    return quad(lambda s: s**m*g1(s)**2*g1p(s)**2, [0,mpf('0.5'),1,2,4,8,16,40])
def K1(m):
    return 4*factorial(m)/mpf(2)**(m+1)*(eta(m+1)-eta(m))

I = {}
for b in range(1,8):
    for a in range(0,8-b):
        I[(a,b)] = Iab(a,b)
        print(f"I({a},{b})[w{a+b+1}] = {nstr(I[(a,b)],40)}")
K2 = {m: K2int(m) for m in [4,5,6]}
K3 = {m: K3int(m) for m in [4,6]}
for m,v in K2.items(): print(f"K2({m}) = {nstr(v,40)}")
for m,v in K3.items(): print(f"K3({m}) = {nstr(v,40)}")
pickle.dump({'I':{str(k):nstr(v,80) for k,v in I.items()},
             'K2':{str(k):nstr(v,80) for k,v in K2.items()},
             'K3':{str(k):nstr(v,80) for k,v in K3.items()}}, open('ckpt_integrals90.pkl','wb'))
