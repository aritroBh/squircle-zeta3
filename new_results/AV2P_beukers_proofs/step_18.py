from mpmath import mp, mpf, zeta, log, pi, nprint
mp.dps = 120
N = 60
zet = [mpf(0)]*(N+3)
for k in range(2, N+3):
    zet[k] = zeta(k)

A = [mpf(0)]*(N+2); A[0] = mpf(1)
for m in range(1, N+1):
    j = m-1
    A[m] = (-1)**(j+1)*zet[j+1]*(mpf(2)**j-1)

def ser_log(A, N):
    L = [mpf(0)]*(N+1)
    Am1 = A[:]; Am1[0] = mpf(0)
    Pr = [mpf(0)]*(N+1); Pr[0] = mpf(1)
    for r in range(1, N+1):
        new = [mpf(0)]*(N+1)
        for i in range(N+1):
            if Pr[i] == 0: continue
            for j in range(N+1-i):
                new[i+j] += Pr[i]*Am1[j]
        Pr = new
        sgn = 1 if r % 2 == 1 else -1
        for n in range(1, N+1):
            L[n] += sgn*Pr[n]/r
    return L

def ser_exp(Lc, N):
    E = [mpf(0)]*(N+1); E[0] = mpf(1)
    for n in range(1, N+1):
        E[n] = sum(k*Lc[k]*E[n-k] for k in range(1, n+1))/n
    return E

L = ser_log(A, N)
c = [mpf(0)]*(N+1)
for k in range(2, N+1):
    c[k] = (-1)**k * zet[k]*(2-mpf(2)**k)/k
logphi = [c[k] + L[k] for k in range(N+1)]
phicoef = ser_exp(logphi, N)
for k in range(0, 12):
    print(k, nprint(phicoef[k], 25))
import json
with open('phicoef.json','w') as f:
    json.dump([str(x) for x in phicoef], f)
print("saved phicoef.json")
