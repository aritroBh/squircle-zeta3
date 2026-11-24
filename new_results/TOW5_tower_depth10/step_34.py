import sympy as sp, time
t0=time.time()
d, m = sp.symbols('d m')
i = sp.I
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
R5 = sp.together(T5)
num5, den5 = sp.fraction(R5)
print("reduced denominator of R_5 = total(5):")
print(sp.factor(den5))
n5e = sp.expand(num5)
print("num degree in d:", sp.degree(n5e, d), " total terms:", len(sp.Poly(n5e, d, m).terms()))
print(f"({round(time.time()-t0,1)}s)")
import pickle
pickle.dump({'num5': str(n5e), 'den5': str(sp.expand(den5))}, open('<project>/research/08_theory/TOW5_tower_depth10/R5_raw.pkl','wb'))
