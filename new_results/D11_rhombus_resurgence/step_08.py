# ============ STEP 1 (fast): exact coefficients via explicit power-series recurrences ============
import sympy as sp

N = 12
eta, s, mu, x, d = sp.symbols('eta s mu x d')
gam = sp.symbols('gam')
z = sp.symbols('z2:z16')
def zet(n): return z[n-2]

def T(poly, var, N):                      # truncate series in var to degree N
    p = sp.expand(poly)
    return sum(p.coeff(var, n)*var**n for n in range(0, N+1))

def series_exp(Lc, N):                    # Lc: list of coeffs (sympy exprs), exp(sum Lc[m] v^m)
    Hc = [sp.Integer(1)]
    for m in range(0, N):
        nxt = sp.expand(sum((j+1)*Lc[j+1]*Hc[m-j] for j in range(0, m+1))/(m+1))
        Hc.append(sp.expand(nxt))
    return Hc

# ---- log G(mu,s):  L(v) with bookkeeping v marking total degree in (mu,s) ----
# Lc[m] = homogeneous part of total degree m of log G(mu,s)
Lc = [sp.Integer(0)]*(N+1)
Lc[1] = -2*gam*mu
for n in range(2, N+1):
    for j in range(0, n+1, 2):
        Lc[n] += 2*(-1)**n * zet(n)/n * sp.binomial(n, j) * mu**(n-j) * s**j
Hc = series_exp([sp.expand(c) for c in Lc], N)
Htot = sum(Hc[m] for m in range(N+1))      # H(mu,s) = exp(logG), total-degree truncated

# ---- c_k (exact rationals) via series sqrt of cosh ----
cc = [sp.Integer(1)]                       # g = sqrt(cosh): coeffs of x^(2k); recurrence g*g=cosh
coshc = [1/sp.factorial(2*k) for k in range(N//2+1)]
for k in range(1, N//2+1):
    cc.append(sp.expand((coshc[k] - sum(cc[i]*cc[k-i] for i in range(1, k)))/2))

# ---- 1/Gamma(2+eta) series ----
logKc = [sp.Integer(0)]*(N+1)
logKc[1] = gam
for n in range(2, N+1):
    logKc[n] = -(-1)**n*zet(n)/n
invG2c = series_exp(logKc, N)
for m in range(0, N):                      # multiply by 1/(1+eta)
    invG2c[m+1] = sp.expand(invG2c[m+1] - (invG2c[m] if m+1 <= N else 0))
# careful: product with geometric series: invG2 *= sum(-eta)^j  -> do properly:
tmp = series_exp(logKc, N)
invG2 = [sp.Integer(0)]*(N+1)
for m in range(N+1):
    invG2[m] = sp.expand(sum(tmp[j]*(-1)**(m-j) for j in range(0, m+1)))

# ---- F(eta) = invG2(eta) * sum_k cc[k] eta^(2k) H_{2k}(eta/2) ----
main = [sp.Integer(0)]*(N+1)
for k in range(0, N//2+1):
    H2k = sp.factorial(2*k) * sp.expand(Htot).coeff(s, 2*k)   # polynomial in mu
    H2k_eta = sp.expand(H2k.subs(mu, eta/2))
    for m in range(0, N+1-2*k):
        main[2*k+m] += cc[k]*sp.expand(H2k_eta).coeff(eta, m)
Fc = [sp.Integer(0)]*(N+1)
for m in range(N+1):
    Fc[m] = sp.expand(sum(invG2[j]*main[m-j] for j in range(0, m+1)))

print("gamma present in any f_n?", any(sp.expand(f).has(gam) for f in Fc[:11]))
for n in range(0, 11):
    print(f"f_{n} =", Fc[n])
