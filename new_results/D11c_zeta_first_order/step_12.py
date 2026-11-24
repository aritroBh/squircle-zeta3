import sympy as sp
from sympy import symbols, Rational, factorial, binomial

NMAX = 12
Z = {n: symbols(f'Z{n}') for n in range(2, NMAX+1)}

def ts_mul(A, B, N=NMAX):
    R = {}
    for i, aa in A.items():
        for j, bb in B.items():
            if i + j <= N: R[i+j] = R.get(i+j, 0) + aa*bb
    return {k: sp.expand(v) for k, v in R.items() if v != 0}
def ts_pow(A, p, N=NMAX):
    R = {0: sp.Integer(1)}
    for n in range(1, N+1):
        acc = 0
        for k in range(1, n+1):
            if k in A and (n-k) in R: acc += (p*k-(n-k))*A[k]*R[n-k]
        v = sp.expand(acc/n)
        if v != 0: R[n] = v
    return R
def ts_exp(A, N=NMAX):
    R = {0: sp.Integer(1)}
    for n in range(1, N+1):
        acc = 0
        for k in range(1, n+1):
            if k in A and (n-k) in R: acc += k*A[k]*R[n-k]
        v = sp.expand(acc/n)
        if v != 0: R[n] = v
    return R
def bs_exp(A, N=NMAX):
    R = {(0,0): sp.Integer(1)}
    for tot in range(1, N+1):
        for i in range(0, tot+1):
            j = tot-i
            if i > 0:
                acc = 0
                for (i1,j1), aa in A.items():
                    if i1 >= 1 and (i-i1, j-j1) in R: acc += i1*aa*R[(i-i1,j-j1)]
                v = sp.expand(acc/i)
            else:
                acc = 0
                for (i1,j1), aa in A.items():
                    if i1 == 0 and j1 >= 1 and (0, j-j1) in R: acc += j1*aa*R[(0,j-j1)]
                v = sp.expand(acc/j) if j > 0 else sp.Integer(1)
            if v != 0: R[(i,j)] = v
    return R

coshm1 = {2*k: Rational(1, factorial(2*k)) for k in range(1, NMAX//2+1)}
sqcosh = ts_pow(coshm1, Rational(1,2))
c = {k: sqcosh.get(2*k, 0) for k in range(NMAX//2+1)}
mA = {n: (-1)**(n+1)*Z[n]/n for n in range(2, NMAX+1)}
K = ts_exp(mA)
Bser = {}
for n in range(2, NMAX+1):
    for r in range(0, n+1, 2):
        Bser[(n-r,r)] = Bser.get((n-r,r), 0) + (-1)**n*Z[n]*2*binomial(n,r)/n
Bser = {k: sp.expand(v) for k, v in Bser.items() if v != 0}
Hser = bs_exp(Bser)
h = {}
for k in range(NMAX//2+1):
    for i in range(0, NMAX+1):
        v = Hser.get((i,2*k), 0)
        if v != 0: h[(k,i)] = sp.expand(factorial(2*k)*v/2**i)
S = {}
for k in range(NMAX//2+1):
    for i in range(0, NMAX+1):
        if (k,i) in h: S[2*k+i] = S.get(2*k+i,0) + c[k]*h[(k,i)]
S = ts_mul(K, S)

lam = {}
for n in range(3, 8, 2):
    for j in range(n, NMAX+1):
        lam[(n,j)] = sp.expand(S.get(j,0)).coeff(Z[n],1)

# Correct L_n: lam_{n,j} = s_{j-n}/n + [eta^j] K(eta) * M_n(eta), M_n = -2 sum_k c_k eta^{2k} sum_{r even} binom(2k,r)(n-1)!/(n-r)! 2^{r-n} eta^{n-r} H^{(2k-2r/2)}(eta/2,0)
def lam_formula(n, j):
    val = S.get(j-n, 0)/n
    M = {}
    for k in range(0, NMAX//2+1):
        for r in range(0, min(2*k,n)+1, 2):
            kk = k - r//2
            fac = -2*c[k]*binomial(2*k,r)*Rational(factorial(n-1),factorial(n-r))*Rational(2)**(r-n)
            for i in range(0, NMAX+1):
                if (kk,i) in h:
                    key = 2*k + (n-r) + i
                    M[key] = M.get(key, 0) + fac*h[(kk,i)]
    KM = ts_mul(K, M)
    return sp.expand(val + KM.get(j, 0))

ok = True
for n in [3,5,7]:
    for j in range(n, NMAX+1):
        if sp.simplify(lam_formula(n,j) - lam[(n,j)]) != 0:
            ok = False; print("MISMATCH", n, j)
print("corrected L_n(eta) formula verified:", ok)

# Print lambda table for reference
for n in [3,5]:
    print(f"lam_{{{n},j}}:", {j: lam[(n,j)] for j in range(n, 11)})
