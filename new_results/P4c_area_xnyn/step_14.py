# step_12.py : verify C = C0 + Delta1 + Delta2 numerically, and the Tauberian lemma for T_k
from mpmath import mp, mpf, sqrt, log, nstr, pi, zeta, gamma, sinh, cosh, exp, tanh, inf, quad
mp.dps = 60
L2 = log(2); z2 = zeta(2)

def fW(w):
    w = mpf(w); return 1 + w - sqrt(1+w*w)
def gW(w):
    w = mpf(w); return fW(w)/w

def C_direct(n):
    n = mpf(n); e = 1/n; a = 1 - e/2
    def inner(sig):
        return quad(lambda d: fW(n*exp(-a*sig)*sqrt(2*cosh(e*d))), [0, sig/2, sig])
    L = log(n)
    return quad(lambda sg: exp(-e*sg)*inner(sg), [0, L/2, L, L+30])

def C0(n):
    # ∫_0^N log(N/w) g(w) dw, N = n sqrt 2
    N = n*sqrt(2)
    return quad(lambda w: log(N/w)*gW(w), [0, 1, N])

def Cprime(n):
    n = mpf(n); e = 1/n
    return quad(lambda sg: exp(-e*sg)*sg*fW(n*sqrt(2)*exp(-sg)), [0, log(n), 2*log(n), 5*log(n)])

def C0prime(n):  # C0(n) but as ∫ σ f dσ (no e^{-eσ})
    n = mpf(n)
    return quad(lambda sg: sg*fW(n*sqrt(2)*exp(-sg)), [0, log(n), 2*log(n)])

for n in [mpf(100), mpf(1000)]:
    C = C_direct(n); C0v = C0(n); Cp = Cprime(n); C0p = C0prime(n)
    D1 = C - Cp; D2 = Cp - C0p
    print("n =", n)
    print("  C direct      =", nstr(C, 30))
    print("  C0 (w-int)    =", nstr(C0v, 30), "  C0 (sig-int) =", nstr(C0p, 30), " diff:", nstr(C0v-C0p, 6))
    print("  C' (e^{-esig}) =", nstr(Cp, 30))
    print("  Delta1 = C-C' =", nstr(D1, 20), "  Delta2 = C'-C0 =", nstr(D2, 20))
    print("  C0+D1+D2 - C  =", nstr(C0v+D1+D2-C, 8))

# Tauberian lemma test: T_k(N) = ∫_0^N log^k(N/w) g dw  vs  L^{k+1}/(k+1) + Σ k!/(k-i)! m_i L^{k-i}
# M_g(z) = E_g/(z(1-z)), E_g = Gamma((1+z)/2) Gamma(1-z/2)/sqrt(pi)
def Eg(z):
    z = mpf(z)
    return gamma((1+z)/2)*gamma(1-z/2)/sqrt(pi)
h = mpf('1e-20')
# series of Eg via high-order finite differences -> use sympy instead for exact e_k
import sympy as sp
zz = sp.symbols('zz')
Eg_series = sp.series(sp.gamma((1+zz)/2)*sp.gamma(1-zz/2)/sp.sqrt(sp.pi), zz, 0, 6).removeO().expand()
print("E_g series:", Eg_series)
eks = [mpf(str(Eg_series.coeff(zz, k))) for k in range(6)]
m = [sum(eks[:k+1]) for k in range(6)]   # m_k = sum_{i<=k} e_i
print("e_k:", [nstr(x,12) for x in eks])
print("m_k:", [nstr(x,12) for x in m])
def T_k(N, k):
    N = mpf(N)
    return quad(lambda w: log(N/w)**k * gW(w), [0, 1, N])
def T_k_pred(N, k):
    L = log(mpf(N))
    s = L**(k+1)/mpf(k+1)
    import math
    for i in range(k+1):
        s += mpf(math.factorial(k))/mpf(math.factorial(k-i)) * m[i] * L**(k-i)
    return s
for N in [mpf(1000)]:
    for k in [0,1,2,3]:
        print("N=1000 k=%d: T=%s pred=%s diff=%s" % (k, nstr(T_k(N,k),20), nstr(T_k_pred(N,k),20), nstr(T_k(N,k)-T_k_pred(N,k),8)))
