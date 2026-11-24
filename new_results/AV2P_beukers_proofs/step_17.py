from mpmath import mp, mpf, zeta, euler, log, exp, pi, nprint, quad, inf
mp.dps = 120

N = 60
# series for log phi:
# logphi(w) = sum_{k>=2} (-1)^k zeta(k)(2-2^k) w^k/k + log(1 + sum_{m>=2} a_{m-1} w^m)
# with a_j = (-1)^{j+1} zeta(j+1)(2^j-1),  from Delta(w) = psi(2w+1)-psi(w) = 1/w + sum a_j w^j
zet = [zeta(k) for k in range(0, N+3)]
A = [mpf(0)]*(N+2)   # coefficients of w^m in (1 + sum a_{m-1} w^m): A[0]=1, A[m]=a_{m-1}
A[0] = mpf(1)
for m in range(1, N+1):
    j = m-1
    A[m] = (-1)**(j+1)*zet[j+1]*(2**j-1)
# log(A) series: L = log(A), L' = A'/A -> compute via series division
# use formal log: log(A) = sum_{r>=1} (-1)^{r+1} (A-1)^r / r
B = A[:]  # A as list of coeffs up to N
# compute log via iterative: L_n = L_{n-1} + ...
import itertools
def ser_log(A, N):
    # A[0]=1
    L = [mpf(0)]*(N+1)
    P = [mpf(0)]*(N+1)  # powers of (A-1)
    Am1 = A[:]; Am1[0] = mpf(0)
    # P_r = (A-1)^r
    Pr = [mpf(0)]*(N+1); Pr[0] = mpf(1)
    for r in range(1, N+1):
        # Pr = Pr * Am1
        new = [mpf(0)]*(N+1)
        for i in range(N+1):
            if Pr[i] == 0: continue
            for j in range(N+1-i):
                new[i+j] += Pr[i]*Am1[j]
        Pr = new
        for n in range(1, N+1):
            L[n] += (-1)**(r+1)*Pr[n]/r
    return L
L = ser_log(B, N)
# add the gamma part
c = [mpf(0)]*(N+1)
for k in range(2, N+1):
    c[k] = (-1)**k * zet[k]*(2-2**k)/k
logphi = [c[k] + L[k] for k in range(N+1)]
# exponentiate
def ser_exp(L, N):
    E = [mpf(0)]*(N+1); E[0] = mpf(1)
    for n in range(1, N+1):
        E[n] = sum(k*L[k]*E[n-k] for k in range(1, n+1))/n
    return E
phicoef = ser_exp(logphi, N)
for k in range(0, 12):
    print(k, nprint(phicoef[k], 25))
import json
with open('phicoef.json','w') as f:
    json.dump([str(x) for x in phicoef], f)
print("saved")
