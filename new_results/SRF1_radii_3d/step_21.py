
import mpmath as mp, json, time
from sympy import partitions as spart
from math import comb, factorial

t0 = time.time()
n, K = 3, 60
mp.mp.dps = 220
zeta = {k: mp.zeta(k) for k in range(2, K+2)}

def ser_exp(L, K):
    E = [mp.mpf(0)]*(K+1); E[0] = mp.mpf(1)
    for k in range(1, K+1):
        s = mp.mpf(0)
        for j in range(1, k+1):
            s += j*L[j]*E[k-j]
        E[k] = s/k
    return E
def ser_mul(A, B, K):
    C = [mp.mpf(0)]*(K+1)
    for i,a in enumerate(A):
        if a==0: continue
        for j in range(K+1-i):
            if B[j]!=0: C[i+j] += a*B[j]
    return C
def ser_log(A, K):
    L = [mp.mpf(0)]*(K+1)
    for nn in range(1, K+1):
        s = nn*A[nn]
        for k in range(1, nn):
            s -= k*L[k]*A[nn-k]
        L[nn] = s/nn
    return L

# Lambda_n
Lam = [mp.mpf(0)]*(K+1)
for k in range(2, K+1):
    Lam[k] = -zeta[k]*mp.mpf(n-1)**k*(1-mp.mpf(n)**(1-k))/k
eLam = ser_exp(Lam, K)

# composition sums CS[s], s=0..K
CS = [None]*(K+1)
CS[0] = [mp.mpf(1)]+[mp.mpf(0)]*K
cm = mp.mpf(-(n-1))/n
for s in range(1, K+1):
    acc = [mp.mpf(0)]*(K+1)
    for p in spart(s):
        parts, mults = [], []
        for part, m in sorted(p.items()):
            parts += [part]*m; mults.append(m)
        r = len(parts)
        if r > n: continue
        coef = mp.mpf(factorial(s))
        for t in range(r): coef *= (n-t)
        for m in mults: coef /= factorial(m)
        for l in parts: coef /= factorial(l)
        cvals = [mp.mpf(2*l*n-2*s-(n-1))/n for l in parts]
        c0 = mp.mpf(-2*s-(n-1))/n
        lm = [mp.mpf(0)]*(K+1)
        for k in range(2, K+1):
            br = sum(c**k for c in cvals) + (n-r)*c0**k - n*cm**k
            if br: lm[k] = ((-1)**k)*br/k
        M = ser_exp(lm, K)
        for i in range(K+1): acc[i] += coef*M[i]
    CS[s] = acc

# B_n
Bser = [mp.mpf(1)]+[mp.mpf(0)]*K
hb = mp.mpf(1)  # binom(1/2, j) updated
for j in range(1, K+1):
    hb = hb*(mp.mpf(1,2)-(j-1))/j
    EX = [mp.mpf(0)]*(K+1)
    for s in range(0, j+1):
        fac = mp.mpf((-1)**(j-s))*comb(j,s)*mp.mpf(n)**(j-s)/mp.mpf(n)**j
        for i in range(K+1): EX[i] += fac*CS[s][i]
    for i in range(K+1): Bser[i] += hb*EX[i]

Sser = ser_mul(eLam, Bser, K)          # S_n/S_n(1)/Rat_n in tau
W = ser_log(Sser, K)                    # w_k
Rat = ser_mul([mp.mpf(1),-mp.mpf(1)]+[mp.mpf(0)]*(K-1), [mp.mpf((2)**r) for r in range(K+1)], K)  # (1-tau)/(1-2tau)
Full = ser_mul(Rat, Sser, K)            # S/(4sqrt3) in tau
# d-series: tau = d/(1+d)
etad = [mp.mpf(0)]+[mp.mpf((-1)**(m-1)) for m in range(1, K+1)]
Dser = [mp.mpf(0)]*(K+1); pw = [mp.mpf(1)]+[mp.mpf(0)]*K
for e in range(K+1):
    if e>0: pw = ser_mul(pw, etad, K)
    for r in range(K+1): Dser[r] += Full[e]*pw[r]
json.dump({'tau':[mp.nstr(x,180) for x in Full],'d':[mp.nstr(x,180) for x in Dser],'w':[mp.nstr(x,180) for x in W]},
          open('SRF1_n3_K60.json','w'))
print(f"engine done in {time.time()-t0:.0f}s")
# validate vs exact K=12
Wexact = json.load(open('w_exact_n3_K12.json'))
worst = mp.mpf(0)
for k in range(2,13):
    vex = mp.mpf(0)
    for m,c in Wexact[str(k)].items():
        t = mp.mpf(sp_Rational := __import__('sympy').Rational(c))
        for i in map(int, m.split()): t *= zeta[i]
        vex += t
    worst = max(worst, abs(vex-W[k]))
print("worst |W_num - W_exact|, k=2..12:", mp.nstr(worst,3))
# validate vs H1 K=80 file
h1 = json.load(open('<project>/research/10_frontier/H1_octahedron_end/hi_coeffs_n3_K80.json'))
worst2 = max(abs(mp.mpf(h1['eta'][k])-Full[k]) for k in range(0,61))
worst3 = max(abs(mp.mpf(h1['d'][k])-Dser[k]) for k in range(0,61))
print("worst |tau-series vs H1| k<=60:", mp.nstr(worst2,3), "  worst |d-series vs H1|:", mp.nstr(worst3,3))
print("sample w_50, w_60:", mp.nstr(W[50],8), mp.nstr(W[60],8))
