# ================= STEP 1: exact Taylor coefficients of F(eta) and P/sqrt2 =================
import sympy as sp

N = 12                      # we want f_0..f_10 (and a_0..a_10)
eta, s, mu, T, x, d = sp.symbols('eta s mu T x d')
gam = sp.symbols('gam')     # Euler gamma, kept symbolic
z = sp.symbols('z2:z14')    # z(0)=zeta(2), ..., z(11)=zeta(13)
def zet(n): return z[n-2]

def trunc_eta(poly, N):
    return sp.expand(poly).series(eta, 0, N+1).removeO().expand()

# ---- log G(mu,s) = logGamma(1+mu+s)+logGamma(1+mu-s)  (total-degree series in mu,s) ----
# logGamma(1+X) = -gam*X + sum_{n>=2} (-1)^n zeta(n) X^n / n
L = -2*gam*mu
for n in range(2, N+1):
    for j in range(0, n+1, 2):
        L += 2*(-1)**n * zet(n)/n * sp.binomial(n, j) * mu**(n-j) * s**j
L = sp.expand(L)
# H(mu,s) = exp(L), truncated to total degree N (bookkeeping variable T)
LT = sp.expand(L.subs({mu: T*mu, s: T*s}))
HT = sp.exp(LT).series(T, 0, N+1).removeO().expand()
H = sp.expand(HT.subs(T, 1))

# ---- c_k : sqrt(cosh x) = sum c_k x^(2k), exact rationals ----
coshser = sum(x**(2*k)/sp.factorial(2*k) for k in range(0, N//2+2))
gser = sp.series(sp.sqrt(coshser), x, 0, N+1).removeO().expand()
ck = [sp.Rational(gser.coeff(x, 2*k)) for k in range(N//2+1)]
print("c_k =", ck[:8])

# ---- 1/Gamma(2+eta) = exp(gam*eta - sum (-1)^n zeta(n) eta^n/n) / (1+eta) ----
logK = gam*eta - sum((-1)**n * zet(n)/n * eta**n for n in range(2, N+1))
invG2 = trunc_eta(sp.exp(logK).series(eta, 0, N+1).removeO(), N) * trunc_eta(1/(1+eta), N)
invG2 = trunc_eta(invG2, N)

# ---- F(eta) = sum_k c_k eta^(2k) [d_s^(2k) G(mu,s)]_0 / Gamma(2+eta),  mu = eta/2 ----
main = 0
for k in range(0, N//2+1):
    H2k = sp.factorial(2*k) * H.coeff(s, 2*k)          # polynomial in mu (contains gam)
    main += ck[k] * eta**(2*k) * H2k.subs(mu, eta/2)
main = trunc_eta(main, N)
Fser = trunc_eta(invG2 * main, N)
f_exact = [sp.expand(Fser).coeff(eta, n) for n in range(0, 11)]
print("\nEuler gamma left in f_n? ", [sp.expand(f).coeff(gam,1) != 0 for f in f_exact])
for n in range(0, 11):
    print(f"f_{n} =", f_exact[n])
