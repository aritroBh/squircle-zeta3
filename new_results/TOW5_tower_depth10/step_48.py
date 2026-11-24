import sympy as sp, pickle, time
t0=time.time()
work = '<project>/research/08_theory/TOW5_tower_depth10'
import sys; sys.path.insert(0, work)
d, m = sp.symbols('d m')
i = sp.I
raw = pickle.load(open(work+'/R5_raw.pkl','rb'))
num5 = sp.sympify(raw['num5'], locals={'d':d,'m':m,'I':i})
den5 = sp.sympify(raw['den5'], locals={'d':d,'m':m,'I':i})
g = sp.gcd(num5, den5)
N5 = sp.expand(sp.div(num5, g)[0])
Pi5 = sp.expand(den5 / g / 3)   # should be Pi_5
print("den5/gcd/3 == Pi_5:", sp.expand(Pi5 - sp.prod(2*d-2*j-1 for j in range(1,6)))==0)
N5c = sp.Poly(N5, d, m)
print("N_5: degree in d:", sp.degree(N5, d), " terms:", len(N5c.terms()))
# check if N_5 divisible by 3 in Z[i][d,m]
coeffs = [c for _, c in N5c.terms()]
from sympy import gcd as igcd
allcoeffs = []
for c in coeffs:
    allcoeffs += [sp.re(c), sp.im(c)]
g3 = sp.Integer(0)
for c in allcoeffs:
    g3 = sp.gcd(g3, sp.Integer(c))
print("gcd of all N_5 coefficients (in Z[i]):", g3)
# R_5 = N_5/(3 Pi_5): verify equals total(5) -- recompute total(5) quickly
def ff(x, c):
    r = sp.Integer(1)
    for u in range(c): r *= (x-u)
    return r
def tf(A):
    return (i/2)**(A//2) if A % 2 == 0 else (1+i)/2 * (i/2)**((A-1)//2)
def tf1(A):
    return (-i/2)**(A//2) if A % 2 == 0 else (1-i)/2 * (-i/2)**((A-1)//2)
x = sp.Symbol('x')
P = {0: sp.Integer(1), 1: x*(x-1)/2, 2: x*(x-1)*(x-2)*(3*x-5)/24,
     3: x*(x-1)*(x-2)**2*(x-3)**2/48,
     4: x*(x-1)*(x-2)*(x-3)*(x-4)*(15*x**3-150*x**2+485*x-502)/5760,
     5: x*(x-1)*(x-2)*(x-3)*(x-4)**2*(x-5)**2*(3*x**2-23*x+38)/11520}
def normalized_summand(l, n, r, a, c):
    b = l - a
    jj = d - n; s = m - r
    t = sp.Integer(1)
    for rr in range(0, n):
        R = d - n + 1 + rr
        t *= -2*R/(2*R-3)
    t *= ff(m, r) * ff(d-m, n-r) / ff(d, n)
    t *= P[r].subs(x, m) / ff(m, r)
    t *= P[n-r].subs(x, d-m) / ff(d-m, n-r)
    t *= sp.binomial(l, a) * sp.Rational(2)**n * sp.binomial(n, a-n) * i**(2*n-a) * sp.factorial(a)
    A1 = 2*(r-n) - c; A2 = -2*r - b + c
    t *= sp.binomial(b, c) * ff(2*(jj-s), c) * (-1)**(b-c) * ff(2*s, b-c) * tf(A1) * tf1(A2)
    return sp.expand(t)
def total(l):
    tot = sp.Integer(0)
    for n in range(0, l+1):
        for r in range(0, n+1):
            for a in range(n, l+1):
                for c in range(0, l-a+1):
                    tot += normalized_summand(l, n, r, a, c)
    return sp.expand(tot)
T5 = total(5)
rem = sp.cancel(T5 - N5/(3*sp.prod(2*d-2*j-1 for j in range(1,6))))
print("total(5) - N_5/(3 Pi_5) == 0:", rem==0, f"({round(time.time()-t0,1)}s)")
