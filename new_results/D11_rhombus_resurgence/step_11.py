# ======== STEP 2a: f_n numerically (mpmath truncated series arithmetic, dps=140), n <= 60 ========
from mpmath import mp, mpf, zeta, factorial, binomial, gamma, pi, sqrt
import json
mp.dps = 140
NMAX = 62

# zeta values
Zv = [mpf(0)]*2 + [zeta(n) for n in range(2, NMAX+3)]

# A_j(mu) = [s^j] log G(mu,s) (GAMMA TERMS EXCLUDED; they cancel with 1/Gamma later)
# A_0[m] = 2*(-1)^m zeta(m)/m (m>=2);  A_j[m] = 2*(-1)^(j+m) zeta(j+m)/(j+m) * C(j+m, j), j even >=2, m>=0
A = {}
for j in range(0, NMAX+1, 2):
    arr = [mpf(0)]*(NMAX+1)
    for m in range(0, NMAX+1-j):
        n = j + m
        if n < 2: continue
        arr[m] = 2*((-1)**n)*Zv[n]/n*binomial(n, j)
    A[j] = arr

def smult(a, b, N):            # multiply mu-series
    c = [mpf(0)]*(N+1)
    for i in range(N+1):
        if a[i] == 0: continue
        ai = a[i]
        for j in range(N+1-i):
            if b[j] != 0: c[i+j] += ai*b[j]
    return c

# H = exp(sum_j s^j A_j) in s (even powers only), coeffs = mu-series
Hs = {0: [mpf(1)] + [mpf(0)]*NMAX}
for j in range(2, NMAX+1, 2):
    acc = [mpf(0)]*(NMAX+1)
    for i in range(2, j+1, 2):
        h = Hs[j-i]
        sc = smult(A[i], h, NMAX)
        for m in range(NMAX+1): acc[m] += i*sc[m]
    Hs[j] = [v/j for v in acc]

# c_k exact rationals via recurrence
from fractions import Fraction
ccf = [Fraction(1)]
for k in range(1, NMAX//2+1):
    ccf.append((Fraction(1, factorial(2*k)) - sum(ccf[i]*ccf[k-i] for i in range(1, k)))/2)
cc = [mpf(c.numerator)/mpf(c.denominator) for c in ccf]

# main(eta) = sum_k c_k eta^(2k) H_{2k}(eta/2) * (2k)!
main = [mpf(0)]*(NMAX+1)
for k in range(0, NMAX//2+1):
    hser = Hs[2*k]
    for m in range(0, NMAX+1-2*k):
        main[2*k+m] += cc[k]*hser[m]*factorial(2*k)/2**m

# K(eta) = exp(-sum_{n>=2} (-1)^n zeta(n) eta^n/n)  [gamma cancels against G's -2 gamma mu]
logK = [mpf(0)]*(NMAX+1)
for n in range(2, NMAX+1): logK[n] = -((-1)**n)*Zv[n]/n
Ks = [mpf(1)] + [mpf(0)]*NMAX
for m in range(1, NMAX+1):
    Ks[m] = sum(j*logK[j]*Ks[m-j] for j in range(2, m+1))/m
# 1/(1+eta)
inv1p = [(-1)**m for m in range(NMAX+1)]
# F = K * inv1p * main
def conv(a, b, N):
    c = [mpf(0)]*(N+1)
    for i in range(N+1):
        ai = a[i]
        if ai == 0: continue
        for j in range(N+1-i): c[i+j] += ai*b[j]
    return c
Fnum = conv(conv(Ks, inv1p, NMAX), main, NMAX)

# a_n: P/sqrt2 = 4(1+eta)F; eta(d) = -d/(1+d)
Peta = [4*(Fnum[m] + (Fnum[m-1] if m else 0)) for m in range(NMAX+1)]
etapow = [[mpf(1)]+[mpf(0)]*NMAX]
cur = [mpf(1)]+[mpf(0)]*NMAX
for m in range(1, NMAX+1):
    cur = conv(cur, [mpf(0), mpf(-1), mpf(1)] + [mpf(0)]*(NMAX-2), NMAX)  # multiply by -d+d^2-... ? 
# build eta(d) powers properly: eta = sum_{j>=1} (-1)^j d^j
eta_d = [mpf(0)] + [(-1)**j for j in range(1, NMAX+1)]
etapow = [[mpf(1)]+[mpf(0)]*NMAX]
for m in range(1, NMAX+1):
    etapow.append(conv(etapow[-1], eta_d, NMAX))
Anum = [mpf(0)]*(NMAX+1)
for m in range(NMAX+1):
    for n in range(m, NMAX+1):
        Anum[n] += Peta[m]*etapow[m][n]

json.dump({'f': [mp.nstr(v, 120) for v in Fnum], 'a': [mp.nstr(v, 120) for v in Anum]}, open('step2_num.json','w'))
print("n : f_n (60 digits shown)")
for n in range(0, 13):
    print(f"{n:2d}  {mp.nstr(Fnum[n], 60)}")
print("...")
for n in [20, 30, 40, 50, 60]:
    print(f"{n:2d}  {mp.nstr(Fnum[n], 40)}")
print()
print("a_n numeric:")
for n in range(0, 11):
    print(f"a_{n} = {mp.nstr(Anum[n], 60)}")
