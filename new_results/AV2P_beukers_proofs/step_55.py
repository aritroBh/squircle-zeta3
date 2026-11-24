from mpmath import mp, mpf, zeta, log, pi, nprint, quad, clsin, binomial, sqrt
mp.dps = 50
z3 = zeta(3); L2 = log(2); C13 = clsin(2,pi/3)
def f(x): return 2*mp.asin(sqrt(x)/2)**2
A = quad(lambda x: f(x)/(1-x), [0,1])
B = quad(lambda x: f(x*x)/(1-x), [0,1])
def H(n):
    n = int(n)
    return sum(1/mpf(k) for k in range(1,n+1)) if n>0 else mpf(0)
T = mp.nsum(lambda n: (H(2*n)-H(n))/(mpf(n)**2*binomial(2*n,n)), [1, mp.inf])
print("A =", nprint(A, 30))
print("B =", nprint(B, 30))
print("A-B =", nprint(A-B, 30))
print("T =", nprint(T, 30))
V = -2*quad(lambda u: mp.asin(u/2)**2/(1+u), [0,1])
print("-2 int asin^2/(1+u) =", nprint(V, 30))
L = quad(lambda t: t*log(1+2*mp.sin(t)), [0, pi/6])
print("L =", nprint(L, 30))
cand = z3/2 - (pi/6)*C13 + pi**2*L2/72
print("L - cand:", nprint(L - cand, 5))
# PSLQ L
from mpmath import pslq, catalan, polylog
G = catalan
basis = [1, pi, pi**2, pi**3, z3, L2, L2**2, pi*L2, C13, G, pi*G, pi*C13, C13*L2, log(3), pi*log(3)]
names = ['1','pi','pi2','pi3','z3','log2','log2^2','pilog2','Cl13','G','piG','piCl13','Cl13log2','log3','pilog3']
r = pslq([L]+basis, tol=mpf('1e-45'), maxcoeff=10**12)
print("L pslq:", r)
