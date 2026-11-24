
import mpmath as mp
mp.mp.dps = 60
pi = mp.pi
def cosh(z): return mp.cosh(z)
def gfun(z, q): return mp.sqrt(mp.cosh(z) + q*mp.sinh(z))
def Hfun(z, q): return mp.exp(-z/2)*(gfun(z,q)+gfun(-z,q))
def Kfun(eta, u): return (1+mp.e**(-u))**(-2-eta)

def I_direct(eta, q):
    # I(eta) = int_0^1 [t(1-t)]^{eta/2} g(eta*l) dt, l = log(t/(1-t))
    f = lambda t: (t*(1-t))**(eta/2) * gfun(eta*mp.log(t/(1-t)), q)
    return mp.quad(f, [0, mp.mpf('0.25'), mp.mpf('0.5'), mp.mpf('0.75'), 1])

def I_laplace(eta, q):
    # I(eta) = int_0^inf e^{-u} H(eta*u) K(eta,u) du
    f = lambda u: mp.exp(-u)*Hfun(eta*u, q)*Kfun(eta,u)
    return mp.quad(f, [0, 1, 5, 20, 100])

print("=== Identity I(eta) = int e^{-u} H(eta u) K(eta,u) du ===")
for q, name in [(mp.mpf(0), "square q=0"), (mp.mpf('0.6'), "rhombus q=0.6 (A^2:B^2=4:1)")]:
    for eta in [mp.mpf('0.1'), mp.mpf('0.3'), mp.mpf('-0.1'), mp.mpf('-0.3'), mp.mpf('0.7')]:
        a = I_direct(eta, q); b = I_laplace(eta, q)
        print(f"{name} eta={eta}: I_direct={mp.nstr(a,30)}  I_lap={mp.nstr(b,30)}  diff={mp.nstr(abs(a-b),3)}")

print()
print("=== H_n exact (rational) and H^{[0]}(x) = sum H_m x^m/m! on real axis ===")
from fractions import Fraction
from math import factorial, gcd
N = 400
# cosh z = sum a_k z^{2k}, a_k = 1/(2k)!
# sqrt series s(z): s_0=1, s_k = (a_k - sum_{i=1}^{k-1} s_i s_{k-i})/2
K = N//2+2
a = [Fraction(1, factorial(2*k)) for k in range(K)]
s = [Fraction(0)]*K
s[0] = Fraction(1)
for k in range(1, K):
    acc = a[k] - sum(s[i]*s[k-i] for i in range(1, k))
    s[k] = acc/2
# H(z) = 2 e^{-z/2} s(z): H_n = 2 sum_{2k<=n} s_k (-1/2)^{n-2k}/(n-2k)!
H = [Fraction(0)]*(N+1)
for n in range(N+1):
    acc = Fraction(0)
    for k in range(0, n//2+1):
        acc += s[k] * Fraction((-1)**(n-2*k), 2**(n-2*k)*factorial(n-2*k))
    H[n] = 2*acc
print("H_0..H_8:", [str(v) for v in H[:9]])
print("H_0..H_8 float:", [float(v) for v in H[:9]])

# H^{[0]}(x) = sum H_m x^m/m!  (convergent, entire)
def Hborel0(x):
    tot = mp.mpf(0)
    for m in range(N+1):
        tot += mp.mpf(H[m].numerator)/mp.mpf(H[m].denominator) * x**m/mp.factorial(m)
    return tot
print()
for x in [5, 10, 20, 40, 80]:
    v = Hborel0(mp.mpf(x))
    print(f"H^[0]({x}) = {mp.nstr(v, 20)}   |.|*x^1.5 = {mp.nstr(abs(v)*x**mp.mpf('1.5'), 6)}")
