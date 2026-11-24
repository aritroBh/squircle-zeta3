# ======== STEP 2b: f_n at dps=200, n<=64; save; compare vs dps=140 run ========
from mpmath import mp, mpf, zeta, binomial, pi
from math import factorial
from fractions import Fraction
import json
mp.dps = 200
NMAX = 64
Zv = [mpf(0)]*2 + [zeta(n) for n in range(2, NMAX+3)]
A = {}
for j in range(0, NMAX+1, 2):
    arr = [mpf(0)]*(NMAX+1)
    for m in range(0, NMAX+1-j):
        n = j + m
        if n >= 2: arr[m] = 2*((-1)**n)*Zv[n]/n*binomial(n, j)
    A[j] = arr
def conv(a, b, N):
    c = [mpf(0)]*(N+1)
    for i in range(N+1):
        ai = a[i]
        if ai == 0: continue
        for j in range(N+1-i):
            if b[j] != 0: c[i+j] += ai*b[j]
    return c
H0 = [mpf(1)] + [mpf(0)]*NMAX
for m in range(1, NMAX+1): H0[m] = sum(i*A[0][i]*H0[m-i] for i in range(2, m+1))/m
Hs = {0: H0}
for j in range(2, NMAX+1, 2):
    acc = [mpf(0)]*(NMAX+1)
    for i in range(2, j+1, 2):
        sc = conv(A[i], Hs[j-i], NMAX)
        for m in range(NMAX+1): acc[m] += i*sc[m]
    Hs[j] = [v/j for v in acc]
ccf = [Fraction(1)]
for k in range(1, NMAX//2+1):
    ccf.append((Fraction(1, factorial(2*k)) - sum(ccf[i]*ccf[k-i] for i in range(1, k)))/2)
cc = [mpf(v.numerator)/mpf(v.denominator) for v in ccf]
main = [mpf(0)]*(NMAX+1)
for k in range(0, NMAX//2+1):
    for m in range(0, NMAX+1-2*k):
        main[2*k+m] += cc[k]*Hs[2*k][m]*factorial(2*k)/mpf(2)**m
logK = [mpf(0)]*(NMAX+1)
for n in range(2, NMAX+1): logK[n] = -((-1)**n)*Zv[n]/n
Ks = [mpf(1)] + [mpf(0)]*NMAX
for m in range(1, NMAX+1): Ks[m] = sum(j*logK[j]*Ks[m-j] for j in range(2, m+1))/m
inv1p = [mpf((-1)**m) for m in range(NMAX+1)]
Fnum = conv(conv(Ks, inv1p, NMAX), main, NMAX)
json.dump({'f200': [mp.nstr(v, 180) for v in Fnum]}, open('step2_num200.json','w'))
old = json.load(open('step2_num.json'))
print("agreement dps140 vs dps200 at n=60:", mp.nstr(abs(mpf(old['f'][60])-Fnum[60])/abs(Fnum[60]), 3))
print("agreement at n=40:", mp.nstr(abs(mpf(old['f'][40])-Fnum[40])/abs(Fnum[40]), 3))
print("f_60 =", mp.nstr(Fnum[60], 65))
