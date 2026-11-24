# R106: Beukers' integral on the quarter squircle.
from mpmath import mp, mpf, gamma, loggamma, digamma, zeta, pi, log, quad, nsum, inf, nstr, exp, catalan, taylor, fsum, expm1, polylog
mp.dps = 40
def phi(w):
    return exp(3*log(w) + 2*loggamma(w) - loggamma(2*w+1))*(digamma(2*w+1)-digamma(w))
def Bsum(p):
    p = mpf(p)
    return nsum(lambda n: phi(n/p)/n**3, [1, inf])
def Bquad(p):
    p = mpf(p)
    # -1/2 iint log(xy)/(1-xy) over x,y>=0, x^p+y^p<=1, inner y-integral numerically
    def inner(x):
        Y = (1-x**p)**(1/p)
        return quad(lambda y: -log(x*y)/(1-x*y)/2, [0, Y])
    return quad(inner, [0, mpf(1)/2**(1/p), 1])
mp.dps = 25
for p in [1, 2, 3]:
    s = Bsum(p); q = Bquad(p)
    print("p=%s sum=%s |sum-quad|=%s" % (p, nstr(s, 22), nstr(abs(s-q), 3)))
mp.dps = 60
print("B(1) - 2zeta3/3 =", nstr(Bsum(1)-2*zeta(3)/3, 3))
print("B(2) - (zeta3/12 + pi G/3) =", nstr(Bsum(2)-(zeta(3)/12+pi*catalan/3), 3))
print("B(1) =", nstr(Bsum(1), 30), "  B(2) =", nstr(Bsum(2), 30))
# Taylor of phi at 0 (removable): use phi(w) = exp(3log w + 2 loggamma(w) - loggamma(2w+1)) (psi(2w+1)-psi(w))
mp.dps = 60
tc = taylor(phi, mpf('1e-30'), 0)  # value
def phiser(w, K=60):
    # w^3 Gamma(w)^2/Gamma(2w+1) = w*Gamma(1+w)^2/Gamma(1+2w); psi(2w+1)-psi(w) = psi(2w+1)-psi(1+w)+1/w
    pass
# exact Taylor via logGamma series: Gamma(1+w)^2/Gamma(1+2w) = exp(sum_{k>=2} (-1)^k zeta(k)(2-2^k)/k w^k)
from mpmath import mpf as F
K = 8
# series arithmetic
def ser_exp(a):
    r = [F(0)]*K; r[0] = F(1)
    for k in range(1, K): r[k] = fsum(j*a[j]*r[k-j] for j in range(1, k+1))/k
    return r
a = [F(0)]*K
for k in range(2, K): a[k] = (-1)**k*zeta(k)*(2-2**k)/k
G = ser_exp(a)
# psi(1+x) = -gamma + sum_{k>=2} (-1)^k zeta(k) x^{k-1}
# w*(psi(2w+1)-psi(1+w)) = sum_{k>=2} (-1)^k zeta(k) (2^{k-1}-1) w^k
D = [F(0)]*K; D[0] = F(1)
for k in range(2, K): D[k] = (-1)**k*zeta(k)*(2**(k-1)-1)
# phi = G * (w psi-diff + 1) = G*D  (since w^3 Gamma(w)^2/Gamma(2w+1) = w Gamma(1+w)^2/Gamma(1+2w))
ph = [fsum(G[i]*D[k-i] for i in range(k+1)) for k in range(K)]
print("phi Taylor coeffs:", [nstr(x, 12) for x in ph[:6]])
print(" -zeta3 =", nstr(-zeta(3), 12), " 9/4 zeta4 =", nstr(F(9)/4*zeta(4), 12))
# eps^2 coefficient: int_0^inf (phi-1) w^-3 dw  ; near 0 use series (phi-1 = O(w^3))
def phim1(w):
    if w < F('0.05'):
        return fsum(ph[k]*w**k for k in range(1, K)) if False else None
    return phi(w)-1
# longer series for small w
K = 70
a = [F(0)]*K
for k in range(2, K): a[k] = (-1)**k*zeta(k)*(2-2**k)/k
G = ser_exp(a)
D = [F(0)]*K; D[0] = F(1)
for k in range(2, K): D[k] = (-1)**k*zeta(k)*(2**(k-1)-1)
ph = [fsum(G[i]*D[k-i] for i in range(k+1)) for k in range(K)]
def g(w):
    if w < F('0.05'): return fsum(ph[k]*w**(k-3) for k in range(3, K))
    return (phi(w)-1)/w**3
c2 = quad(g, [0, F('0.05'), F(1)/2, 1, 2, 4, 8, 16, 64, inf])
print("int (phi-1)/w^3 =", nstr(c2, 40), "  -pi^2/12 =", nstr(-pi**2/12, 40), " diff", nstr(abs(c2+pi**2/12), 3))
# check expansion against the sum at large p
mp.dps = 40
for p in [100, 200, 400]:
    e = F(1)/p
    B = Bsum(p)
    r = (B - zeta(3) + pi**2/12*e**2 - zeta(3)/2*e**3 + 3*zeta(4)/16*e**4)/e**5
    print("p=%d (B - [zeta3 - pi^2/12 e^2 + zeta3/2 e^3 - 3zeta4/16 e^4])/e^5 = %s" % (p, nstr(r, 10)))
print("predicted e^5 coeff: -(-ph5)*zeta(-2)=0 -> e^5 coeff is 0; e^6: -(-ph6)... check: ph5 =", nstr(ph[5], 10), " zeta(-2)=0; ph6*zeta(-3) =", nstr(ph[6]*zeta(-3), 10))
print("(B - zeta3 + pi^2/12 e^2)/e^3 at p=400:", nstr((Bsum(400)-zeta(3)+pi**2/12/F(400)**2)*400**3, 10), " printed 0.6005 ; zeta3/2 =", nstr(zeta(3)/2, 8))
