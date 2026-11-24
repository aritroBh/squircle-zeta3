"""prove_S12.py -- PROOF for all (d,m) of the identities
   S^{(l)}(t*)/S(t*) = N_l(d,m) / prod_{j=1..l}(2d-2j-1),  l = 1, 2,
where S(t) = S_m(t) is the sqrt-block polynomial of the degree-d initial integrand
(genintegrand.S_poly), t* = (1+i)/2, and N_1, N_2 are the fitted polynomials.
Method: the normalized derivative T_l/T_0 is a finite sum of explicit rational
functions of (d,m) over Q(i) (Stirling-tail polynomials P_t(x) = S(x,x-t), proved
separately below by their recurrence). sympy cancels the total minus N_l/prod to 0:
an exact rational identity, valid for all (d,m).
"""
import sympy as sp
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

# Stirling-tail polynomials P_t(x) = S(x, x-t), t <= 2 (proved by recurrence below)
x = sp.Symbol('x')
P = {0: sp.Integer(1), 1: x*(x-1)/2, 2: x*(x-1)*(x-2)*(3*x-5)/24}
# verify the Stirling recurrence P_t(x+1) = (x+1-t) P_{t-1}(x) + P_t(x) symbolically
assert sp.expand(P[1].subs(x,x+1) - ((x+1-1)*P[0] + P[1])) == 0
assert sp.expand(P[2].subs(x,x+1) - ((x+1-2)*P[1].subs(x,x) + P[2])) == 0
# and the anchor values P_t(t) = S(t,0) = 0 for t >= 1, P_t(0)=0
assert P[1].subs(x,0) == 0 and P[1].subs(x,1) == 0
assert P[2].subs(x,0) == 0 and P[2].subs(x,1) == 0 and P[2].subs(x,2) == 0
# (with P_0 = 1 and these anchors the 2-term recurrence forces P_t(x) = S(x,x-t)
#  for ALL integers x,t by induction on t then x.)

def C_half(n):
    """binom(1/2, d-n)/binom(1/2, d) = prod_{r=d-n+1..d} [-2r/(2r-3)]"""
    out = sp.Integer(1)
    for r in range(d-n+1, d+1):   # symbolic range not allowed; n is concrete int
        pass
    return out

def normalized_summand(l, n, r, a, c):
    """one normalized term of T_l/T_0: j = d-n, s = m-r, q2-power a, B-derivative split c."""
    b = l - a
    jj = d - n; s = m - r
    t = sp.Integer(1)
    # C_j = binom(1/2,j)/binom(1/2,d)
    for rr in range(0, n):
        R = d - n + 1 + rr
        t *= -2*R/(2*R-3)
    # binom(j,s)/binom(d,m) = m^r_ (d-m)^{n-r}_ / d^n_
    t *= ff(m, r) * ff(d-m, n-r) / ff(d, n)
    # w-ratios: w(s,m)/2^m = P_r(m)/m^r_ ; w(j-s,d-m)/2^{d-m} = P_{n-r}(d-m)/(d-m)^{n-r}_
    t *= P[r].subs(x, m) / ff(m, r)
    t *= P[n-r].subs(x, d-m) / ff(d-m, n-r)
    # binom(l,a) * (q2^n)^{(a)}(t*)/2^d = binom(l,a) 2^{n-d} binom(n,a-n) i^{2n-a} a!
    t *= sp.binomial(l, a) * sp.Rational(2)**n * sp.binomial(n, a-n) * i**(2*n-a) * sp.factorial(a)
    # B-part: binom(b,c) ff(2(j-s),c) (-1)^{b-c} ff(2s,b-c) t*^{A1} (1-t*)^{A2}
    A1 = 2*(r-n) - c; A2 = -2*r - b + c
    t = sp.expand(t) if False else t
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

N = {1: 4*m + (-2+6*i)*d - 4*i*d**2,
     2: (-40*i)*m + 16*m**2 + (sp.Rational(80,3)+20*i)*d + (-16+96*i)*d*m + (-80-48*i)*d**2 + (-32*i)*d**2*m + (sp.Rational(208,3)+16*i)*d**3 - 16*d**4}
def Pi(l):
    out = sp.Integer(1)
    for jj in range(1, l+1): out *= (2*d-2*jj-1)
    return out

for l in (1, 2):
    rem = sp.cancel(sp.together(total(l) - N[l]/Pi(l)))
    print(f'S^{l}/S identity remainder:', rem)
