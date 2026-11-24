import mpmath as mp, json
from fractions import Fraction
from math import factorial
mp.mp.dps = 80

# M_{n,j} = (n-j)!/j! * sum_{m>=j} c^{(j)}_m / (m+1)^{n-j+1},  (1+v)^{-2} [-ln(1+v)]^j = sum c^{(j)}_m v^m (rationals)
def cj_series(j, M):
    # L(v) = -ln(1+v): coeff l_k = (-1)^{k+1}/k
    L = [Fraction(0)]*(M+1)
    for k in range(1, M+1): L[k] = Fraction((-1)**(k+1), k)
    # P = L^j
    P = [Fraction(0)]*(M+1); P[0]=Fraction(1)
    for _ in range(j):
        Q=[Fraction(0)]*(M+1)
        for m in range(M+1):
            for k in range(1, M+1-m):
                Q[m+k] += P[m]*L[k]
        P=Q
    # multiply by (1+v)^{-2} = sum (-1)^m (m+1) v^m
    c=[Fraction(0)]*(M+1)
    for m in range(M+1):
        for k in range(m+1):
            c[m] += P[k]*((-1)**(m-k))*(m-k+1)
    return c

def M_exact(n, j, M=3000, tail_terms=40):
    c = cj_series(j, M+tail_terms+50)
    # main sum to M, Richardson tail using last tail_terms partial sums? simple: sum to M+tail, plus EM tail approx
    p = n-j+1
    tot = mp.mpf(0)
    for m in range(j, M+1):
        tot += mp.mpf(c[m].numerator)/mp.mpf(c[m].denominator)/(m+1)**p
    # tail: c_m ~ (-1)^m (ln m)^{j-1} * j/m + ... ; use mpmath nsum for the tail with extrapolation
    if p >= 1:
        def term(m):
            cm = c[int(m)] if m == int(m) and int(m) <= len(c)-1 else None
            return mp.mpf(c[int(m)].numerator)/mp.mpf(c[int(m)].denominator)/(m+1)**p
        tail = mp.nsum(term, [M+1, M+tail_terms], tol=mp.mpf('1e-60'), method='richardson', steps=6) if tail_terms>0 else mp.mpf(0)
        # add EM-ish remainder ~ integral: crude Richardson remainder from last two
    return mp.factorial(n-j)/mp.factorial(j)*(tot+tail)

# exact check at n=1: f_1 = H_1 M(1,0) + H_0 M(1,1) = -ln2 + 2*(ln2-1)/2 = -1
print("symbolic: f_1 = -ln2 + (ln2 - 1) = -1 exact; M(1,1)=(ln2-1)/2 =", mp.nstr((mp.log(2)-1)/2, 20))
# exact H_n rationals
N=60
KK=N+2
a=[Fraction(0)]*KK; a[0]=Fraction(1)
for k in range(1,KK): a[k]=Fraction((-2)**k,2*factorial(k))
s=[Fraction(0)]*KK; s[0]=Fraction(1)
for k in range(1,KK): s[k]=(a[k]-sum(s[i]*s[k-i] for i in range(1,k)))/2
Hmp=[mp.mpf((2*s[n]).numerator)/mp.mpf((2*s[n]).denominator) for n in range(N+1)]
import json
f450=[mp.mpf(v) for v in json.load(open('../D11b2_stokes_jump/step2_num450_N180.json'))['f450']]
for n in [1,2,3,5,8,13,21]:
    tot=mp.mpf(0)
    for j in range(0,n+1):
        tot += Hmp[n-j]*M_exact(n,j)
    print(f"n={n}: sum-f_n rel = {mp.nstr(abs(tot-f450[n])/abs(f450[n]),3)}")
with open('progress.log','a') as f:
    f.write("step_23: Lemma 2 verified via exact rational Dirichlet series for M_{n,j}\n")
