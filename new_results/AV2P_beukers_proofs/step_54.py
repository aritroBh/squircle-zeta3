from mpmath import mp, mpf, zeta, log, pi, nprint, quad, clsin, binomial
mp.dps = 50
z3 = zeta(3); L2 = log(2); C13 = clsin(2,pi/3)
def f(x): return 2*mp.asin(sqrt(x)/2)**2
A = quad(lambda x: f(x)/(1-x), [0,1])
B = quad(lambda x: f(x*x)/(1-x), [0,1])
T = mp.nsum(lambda n: (mp.nsum(lambda j: 1/mpf(n+j), [1,n]))/(mpf(n)**2*binomial(2*n,n)), [1, mp.inf])
print("A =", nprint(A, 30))
print("B =", nprint(B, 30))
print("A-B =", nprint(A-B, 30))
print("T =", nprint(T, 30))
# the -2 int form:
V = -2*quad(lambda u: mp.asin(u/2)**2/(1+u), [0,1])
print("-2 int asin^2/(1+u) =", nprint(V, 30))
# so what is the right closed form of that integral?
# IBP on V: V = -2[asin^2(u/2) log(1+u)]_0^1 + 4 int log(1+u) asin(u/2)/sqrt(4-u^2) du
V2 = -2*(mp.asin(mpf('0.5'))**2)*log(2) + 4*quad(lambda u: log(1+u)*mp.asin(u/2)/sqrt(4-u**2), [0,1])
print("V2 (IBP) =", nprint(V2, 30))
# L = int_0^{pi/6} th log(1+2 sin th) dth: value and PSLQ
L = quad(lambda t: t*log(1+2*mp.sin(t)), [0, pi/6])
print("L =", nprint(L, 30))
cand = z3/2 - (pi/6)*C13 + pi**2*L2/72
print("L - cand:", nprint(L - cand, 5))
print("cand =", nprint(cand, 30))
