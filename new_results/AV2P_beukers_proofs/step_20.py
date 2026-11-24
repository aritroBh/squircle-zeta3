from mpmath import mp, mpf, zeta, log, pi, nprint, quad, inf, exp, loggamma, digamma
mp.dps = 130
N = 250
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

a = mpf('0.2')
# small-w part: sum_{k>=3} phi_k a^{k-2}/(k-2)
small = sum(phicoef[k]*a**(k-2)/(k-2) for k in range(3, N+1))
print("small part:", nprint(small, 40))
def phi_direct(w):
    w = mpf(w)
    return w**3 * exp(2*loggamma(w) - loggamma(2*w+1)) * (digamma(2*w+1) - digamma(w))
ints = [(a, 1), (1, 4), (4, 16), (16, 64), (64, 200), (200, inf)]
tot = small
for aa,b in ints:
    tot += quad(lambda w: (phi_direct(w)-1)/w**3, [aa, b])
C2 = tot
print("C2 =", nprint(C2, 110))
print("C2 + pi^2/12 =", nprint(C2 + pi**2/12, 6))
import json
with open('phicoef.json','w') as f:
    json.dump([str(x) for x in phicoef], f)
with open('progress.log','a') as fh:
    fh.write("step_18: phi Taylor series N=250 computed; C2 integral to ~110 digits matches -pi^2/12\n")
