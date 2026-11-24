# step_17.py : dissect C(n=1e6) : C0, C', Delta1, Delta2, poly series, tail series
from mpmath import mp, mpf, sqrt, log, nstr, pi, zeta, inf, quad, exp, cosh
import math
mp.dps = 70
L2 = log(2); z2 = pi**2/6; z3 = zeta(3); z4 = pi**4/90

def fW(w):
    w = mpf(w); return 1 + w - sqrt(1+w*w)
def gW(w): return fW(w)/w
def fp(w):
    w = mpf(w); return 1 - w/sqrt(1+w*w)

n = mpf(10)**6; e = 1/n; N = n*sqrt(2)
# C0 = T_1(g)(N) exactly
C0 = quad(lambda w: log(N/w)*gW(w), [0, 1, N])
# C' = ∫ e^{-eσ} σ f(W0) dσ
Cp = quad(lambda sg: exp(-e*sg)*sg*fW(n*sqrt(2)*exp(-sg)), [0, log(n), 2*log(n), 4*log(n), 12*log(n)])
D2 = Cp - C0
# Delta1 direct
def D1_direct():
    def inner(sg):
        W0 = n*sqrt(2)*exp(-sg)
        return quad(lambda d: fW(W0*exp(e*sg/2)*sqrt(cosh(e*d))) - fW(W0), [0, sg/2, sg])
    L = log(n)
    return quad(lambda sg: exp(-e*sg)*inner(sg), [0, L/2, L, 2*L, 5*L])
D1 = D1_direct()
C_actual = C0 + D1 + D2
C_fromE = mpf('105.9965931431501241549511127430218278965')  # from E_values_new (n=1e6)
print("C0 exact      =", nstr(C0, 40))
print("C'            =", nstr(Cp, 40))
print("D2 = C'-C0    =", nstr(D2, 30))
print("D1            =", nstr(D1, 30))
print("C = C0+D1+D2  =", nstr(C_actual, 40))
print("C from E      =", nstr(C_fromE, 40))
print("diff          =", nstr(C_actual - C_fromE, 12))

# poly_k(g) formulas
e1=-L2; e2=L2**2/2+z2/2; e3=-z3/4-L2*z2/2-L2**3/6; e4=z4/4+L2*z3/4+L2**2*z2/4+z2**2/8+L2**4/24
eks=[mpf(1),e1,e2,e3,e4]; mt=[sum(eks[:k+1]) for k in range(5)]
def poly(k, N):
    L = log(mpf(N))
    s = L**(k+1)/mpf(k+1)
    for i in range(1, k+1):
        s += mpf(math.factorial(k))/mpf(math.factorial(k+1-i))*mt[i]*L**(k+1-i)
    return s
def r(k, N):  # tail prediction (-1)^k k!/(2N)
    return (-1)**k * mpf(math.factorial(k))/(2*mpf(N))
print()
print("poly1 =", nstr(poly(1,N), 30), "  C0 - poly1 =", nstr(C0-poly(1,N), 12), "  r1 pred =", nstr(r(1,N), 12))
print("poly2 =", nstr(poly(2,N), 30))
# Delta2 poly series
D2poly = sum((-e)**j/mpf(math.factorial(j))*poly(j+1, N) for j in range(1, 12))
D2tail = sum((-e)**j/mpf(math.factorial(j))*r(j+1, N) for j in range(1, 40))
print("D2 actual     =", nstr(D2, 30))
print("D2 polyseries =", nstr(D2poly, 30))
print("D2 tailseries =", nstr(D2tail, 30))
print("D2 - poly     =", nstr(D2 - D2poly, 20))
print("D2 - poly - tail =", nstr(D2 - D2poly - D2tail, 20))
# eT1 check
print("D1 - e*C0 =", nstr(D1 - e*C0, 15))
