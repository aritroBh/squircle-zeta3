import sympy as sp, time
d, m = sp.symbols('d m')
i = sp.I
t0 = time.time()

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
def Pi(l):
    out = sp.Integer(1)
    for jj in range(1, l+1): out *= (2*d-2*jj-1)
    return out
def total(l):
    tot = sp.Integer(0)
    for n in range(0, l+1):
        for r in range(0, n+1):
            for a in range(n, l+1):
                for c in range(0, l-a+1):
                    tot += normalized_summand(l, n, r, a, c)
    return sp.expand(tot)

# Check l=4 works (machinery sanity)
N4 = 13824*i*m - 30784*m**2 - 6912*i*m**3 + 256*m**4 + (-sp.Rational(26576,5)-6912*i)*d + (30784-88320*i)*d*m + (72448+10368*i)*d*m**2 + (-512+6144*i)*d*m**3 + (sp.Rational(60592,3)+44160*i)*d**2 + (-72448+169856*i)*d**2*m + (-54528-9216*i)*d**2*m**2 + (-1024*i)*d**2*m**3 + (-48000-86656*i)*d**3 + (54784-148736*i)*d**3*m + (15872+1536*i)*d**3*m**2 + (sp.Rational(202624,3)+75904*i)*d**4 + (-15872+64768*i)*d**4*m + (-1536)*d**4*m**2 + (-sp.Rational(251904,5)-32640*i)*d**5 + (1536-13312*i)*d**5*m + (sp.Rational(57856,3)+6656*i)*d**6 + (1024*i)*d**6*m + (-3584-512*i)*d**7 + 256*d**8
T4 = total(4)
rem4 = sp.together(T4 - N4/Pi(4))
print("l=4 check: numerator == 0:", sp.expand(sp.numer(rem4))==0, f"({round(time.time()-t0,1)}s)")

# Now: at integer points, compare total(5) vs direct S^(5)/S, and check (S5/S)*Pi5 polynomial
from math import comb, factorial
def stir2(n, k):
    if n == k == 0: return 1
    if n <= 0 or k <= 0 or k > n: return 0
    return sum((-1)**(k-j) * comb(k, j) * j**n for j in range(k+1)) // factorial(k)
def S_poly(dd, pp):
    qq = dd - pp
    t = sp.Symbol('t')
    W = t**2 + (1-t)**2
    out = sp.Integer(0)
    for j in range(1, dd+1):
        Cj = sp.Rational(sp.binomial(sp.Rational(1,2), j) * factorial(j)) * 2**dd
        slo = max(0, j-qq); shi = min(j, pp)
        for s in range(slo, shi+1):
            Sps = stir2(pp, s); Sqj = stir2(qq, j-s)
            if Sps == 0 or Sqj == 0: continue
            out += Cj * Sps * Sqj * t**(2*(j-s)) * (1-t)**(2*s) * W**(dd-j)
    return sp.expand(out) / (factorial(pp)*factorial(qq))
tstar = (1+i)/2
T5 = total(5)
match = 0; mism = 0
for dd in range(6, 11):
    for pp in range((dd+1)//2, dd+1):
        Sp = S_poly(dd, pp)
        direct = sp.simplify(sp.diff(Sp,'t',5).subs('t',tstar) / Sp.subs('t',tstar))
        viam = sp.simplify(T5.subs({d:dd,m:pp}))
        if sp.simplify(direct-viam)==0: match += 1
        else: mism += 1; print("MISMATCH total(5) vs direct at", dd, pp)
print(f"total(5) == direct S5/S at integer points: match={match} mism={mism} ({round(time.time()-t0,1)}s)")
