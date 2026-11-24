"""rhombus_exact_coefficients.py -- exact Taylor coefficients of the squircle perimeter at the diamond p = 1.

Proposition 2 of "Investigating zeta(3)". With eta = 1/p - 1 and mu = eta/2,
    P/(4 sqrt2 (1+eta)) = F(eta) = sum_k c_k eta^(2k) d_s^(2k)[Gamma(1+mu+s) Gamma(1+mu-s)]_(s=0) / Gamma(1+eta),
where sqrt(cosh z) = sum c_k z^(2k). Expanding log Gamma(1+x) = -gamma x + sum (-1)^n zeta(n) x^n / n
symbolically (Z2, Z3, ... stand for zeta(2), zeta(3), ...), this script prints f_0..f_10 and the
coefficients a_0..a_10 of P/sqrt2 = sum a_n d^n, d = p - 1, with even zeta values reduced to powers of pi.
It also prints whether the Euler-Mascheroni constant survives in any coefficient (it does not).
Run: python rhombus_exact_coefficients.py   (a few seconds)
"""
import sympy as sp

N = 12
eta, s, mu, x, d = sp.symbols('eta s mu x d')
gam = sp.symbols('gam')
Z = sp.symbols(' '.join(f'Z{n}' for n in range(2, 15)))
def zet(n): return Z[n-2]
def T(poly, var, N):
    p = sp.expand(poly)
    return sum(p.coeff(var, n)*var**n for n in range(0, N+1))
def series_exp(Lc, N):
    Hc = [sp.Integer(1)]
    for m in range(0, N):
        Hc.append(sp.expand(sum((j+1)*Lc[j+1]*Hc[m-j] for j in range(0, m+1))/(m+1)))
    return Hc

# log G(mu,s) homogeneous parts (total degree in mu,s)
Lc = [sp.Integer(0)]*(N+1); Lc[1] = -2*gam*mu
for n in range(2, N+1):
    for j in range(0, n+1, 2):
        Lc[n] += 2*(-1)**n * zet(n)/n * sp.binomial(n, j) * mu**(n-j) * s**j
Hc = series_exp([sp.expand(c) for c in Lc], N)
Htot = sum(Hc)
# c_k via g*g = cosh recurrence
cc = [sp.Integer(1)]; coshc = [1/sp.factorial(2*k) for k in range(N//2+1)]
for k in range(1, N//2+1):
    cc.append(sp.expand((coshc[k] - sum(cc[i]*cc[k-i] for i in range(1, k)))/2))
# 1/Gamma(2+eta)
logKc = [sp.Integer(0)]*(N+1); logKc[1] = gam
for n in range(2, N+1): logKc[n] = -(-1)**n*zet(n)/n
tmp = series_exp(logKc, N)
invG2 = [sp.expand(sum(tmp[j]*(-1)**(m-j) for j in range(0, m+1))) for m in range(N+1)]
# F(eta)
main = [sp.Integer(0)]*(N+1)
for k in range(0, N//2+1):
    H2k_eta = sp.expand((sp.factorial(2*k)*sp.expand(Htot).coeff(s, 2*k)).subs(mu, eta/2))
    for m in range(0, N+1-2*k):
        main[2*k+m] += cc[k]*sp.expand(H2k_eta).coeff(eta, m)
Fc = [sp.expand(sum(invG2[j]*main[m-j] for j in range(0, m+1))) for m in range(N+1)]
# a_n in d
eta_ser = sum((-d)**j for j in range(1, 13))
eta_pow = [sp.Integer(1)]
for m in range(1, 12): eta_pow.append(T(eta_pow[-1]*eta_ser, d, 11))
Ac_eta = [sp.expand(4*(Fc[m] + (Fc[m-1] if m >= 1 else 0))) for m in range(0, 11)]
Ac = [sp.expand(sum(Ac_eta[m]*eta_pow[m] for m in range(0, 11))).coeff(d, n) for n in range(0, 11)]

pi_red = {zet(2): sp.pi**2/6, zet(4): sp.pi**4/90, zet(6): sp.pi**6/945, zet(8): sp.pi**8/9450, zet(10): sp.pi**10/93555}
out = []
out.append("gamma in f_n? " + str(any(sp.expand(f).has(gam) for f in Fc[:11])))
for n in range(0, 11):
    out.append(f"f_{n} = {sp.expand(Fc[n])}")
for n in range(0, 11):
    out.append(f"a_{n} = {sp.expand(Ac[n].subs(pi_red))}")
txt = "\n".join(out)

print(txt)
