# ======== STEP 2a (fixed): f_n numerically, dps=140 ========
from mpmath import mp, mpf, zeta, binomial, pi
from math import factorial
import json
mp.dps = 140
NMAX = 62
Zv = [mpf(0)]*2 + [zeta(n) for n in range(2, NMAX+3)]
A = {}
for j in range(0, NMAX+1, 2):
    arr = [mpf(0)]*(NMAX+1)
    for m in range(0, NMAX+1-j):
        n = j + m
        if n < 2: continue
        arr[m] = 2*((-1)**n)*Zv[n]/n*binomial(n, j)
    A[j] = arr
def conv(a, b, N):
    c = [mpf(0)]*(N+1)
    for i in range(N+1):
        ai = a[i]
        if ai == 0: continue
        for j in range(N+1-i):
            if b[j] != 0: c[i+j] += ai*b[j]
    return c
Hs = {0: [mpf(1)] + [mpf(0)]*NMAX}
for j in range(2, NMAX+1, 2):
    acc = [mpf(0)]*(NMAX+1)
    for i in range(2, j+1, 2):
        sc = conv(A[i], Hs[j-i], NMAX)
        for m in range(NMAX+1): acc[m] += i*sc[m]
    Hs[j] = [v/j for v in acc]
from fractions import Fraction
ccf = [Fraction(1)]
for k in range(1, NMAX//2+1):
    ccf.append((Fraction(1, factorial(2*k)) - sum(ccf[i]*ccf[k-i] for i in range(1, k)))/2)
cc = [mpf(v.numerator)/mpf(v.denominator) for v in ccf]
main = [mpf(0)]*(NMAX+1)
for k in range(0, NMAX//2+1):
    hser = Hs[2*k]
    for m in range(0, NMAX+1-2*k):
        main[2*k+m] += cc[k]*hser[m]*factorial(2*k)/mpf(2)**m
logK = [mpf(0)]*(NMAX+1)
for n in range(2, NMAX+1): logK[n] = -((-1)**n)*Zv[n]/n
Ks = [mpf(1)] + [mpf(0)]*NMAX
for m in range(1, NMAX+1):
    Ks[m] = sum(j*logK[j]*Ks[m-j] for j in range(2, m+1))/m
inv1p = [mpf((-1)**m) for m in range(NMAX+1)]
Fnum = conv(conv(Ks, inv1p, NMAX), main, NMAX)
Peta = [4*(Fnum[m] + (Fnum[m-1] if m else 0)) for m in range(NMAX+1)]
eta_d = [mpf(0)] + [mpf((-1)**j) for j in range(1, NMAX+1)]
etapow = [[mpf(1)]+[mpf(0)]*NMAX]
for m in range(1, NMAX+1):
    etapow.append(conv(etapow[-1], eta_d, NMAX))
Anum = [mpf(0)]*(NMAX+1)
for m in range(NMAX+1):
    for n in range(m, NMAX+1):
        Anum[n] += Peta[m]*etapow[m][n]
json.dump({'f': [mp.nstr(v, 125) for v in Fnum], 'a': [mp.nstr(v, 125) for v in Anum]}, open('step2_num.json','w'))
print("f_n, n=0..12 (60 digits):")
for n in range(0, 13): print(f"{n:2d}  {mp.nstr(Fnum[n], 60)}")
print("f_n large n (40 digits):")
for n in [20, 30, 40, 50, 60]: print(f"{n:2d}  {mp.nstr(Fnum[n], 40)}")
print("\na_n numeric (60 digits):")
for n in range(0, 11): print(f"a_{n} = {mp.nstr(Anum[n], 60)}")
