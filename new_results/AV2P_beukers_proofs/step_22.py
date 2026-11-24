from mpmath import mp, mpf, zeta, log, pi, nprint, exp, loggamma, digamma, binomial, harmonic, inf, nsum
mp.dps = 220

def phi(w):
    w = mpf(w)
    return w**3 * exp(2*loggamma(w) - loggamma(2*w+1)) * (digamma(2*w+1) - digamma(w))

# B2(2) via Dirichlet series (terms decay 4^{-n/2})
def B2_series(p):
    p = mpf(p); eps = 1/p
    S = mpf(0); n = 1
    nmin = int(500*p)+10
    while True:
        t = phi(n*eps)/mpf(n)**3
        S += t
        if n > nmin and abs(t) < mpf(10)**(-215):
            break
        n += 1
    return S

B2_2 = B2_series(2)
even = B2_series(1)/8   # even part = B2(1)/8 = zeta(3)/12
print("B2(2) =", nprint(B2_2, 60))
print("even part - zeta(3)/12 =", nprint(even - zeta(3)/12, 6))
odd = B2_2 - zeta(3)/12
print("odd part =", nprint(odd, 200))

# odd part via central-binomial series: (pi/8) sum_{k>=0} C(2k,k)/16^k (2ln2+H_k-H_{2k}+1/(2k+1))/(2k+1)
S = mpf(0)
K = 700
Hk = mpf(0)
for k in range(0, K):
    if k > 0:
        Hk += 1/mpf(k)
    H2k = Hk + sum(1/mpf(k+j) for j in range(1, k+1)) if k > 0 else mpf(0)
    term = binomial(2*k, k)/mpf(16)**k * (2*log(2) + Hk - H2k + 1/mpf(2*k+1))/(2*k+1)
    S += term
odd_series = pi/8 * S
print("odd via binomial series =", nprint(odd_series, 60))
print("match:", nprint(odd - odd_series, 6))
with open('odd2.txt','w') as f:
    f.write(str(odd))
with open('progress.log','a') as fh:
    fh.write("step_20: B2(2) even part = B2(1)/8 = zeta(3)/12 PROVED (by Task 1); odd part at 200d via two routes\n")
