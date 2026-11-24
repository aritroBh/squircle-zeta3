import sympy as sp
from sympy import symbols, Rational, factorial, binomial

NMAX = 12
Z = {n: symbols(f'Z{n}') for n in range(2, NMAX+1)}
eta, d = symbols('eta d')

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
            if k in A and (n-k) in R: acc += (p*k - (n-k)) * A[k] * R[n-k]
        val = sp.expand(acc / n)
        if val != 0: R[n] = val
    return R
def ts_exp(A, N=NMAX):
    R = {0: sp.Integer(1)}
    for n in range(1, N+1):
        acc = 0
        for k in range(1, n+1):
            if k in A and (n-k) in R: acc += k * A[k] * R[n-k]
        val = sp.expand(acc / n)
        if val != 0: R[n] = val
    return R

coshm1 = {2*k: Rational(1, factorial(2*k)) for k in range(1, NMAX//2+1)}
sqcosh = ts_pow(coshm1, Rational(1,2))
c = {k: sqcosh.get(2*k, 0) for k in range(NMAX//2+1)}
mA = {n: (-1)**(n+1) * Z[n] / n for n in range(2, NMAX+1)}
K = ts_exp(mA)

def bs_exp(A, N=NMAX):
    R = {(0,0): sp.Integer(1)}
    for tot in range(1, N+1):
        for i in range(0, tot+1):
            j = tot - i
            if i > 0:
                acc = 0
                for (i1,j1), aa in A.items():
                    if i1 >= 1 and (i-i1, j-j1) in R: acc += i1 * aa * R[(i-i1, j-j1)]
                val = sp.expand(acc / i)
            else:
                acc = 0
                for (i1,j1), aa in A.items():
                    if i1 == 0 and j1 >= 1 and (0, j-j1) in R: acc += j1 * aa * R[(0, j-j1)]
                val = sp.expand(acc / j) if j > 0 else sp.Integer(1)
            if val != 0: R[(i,j)] = val
    return R

Bser = {}
for n in range(2, NMAX+1):
    for r in range(0, n+1, 2):
        Bser[(n-r, r)] = Bser.get((n-r, r), 0) + (-1)**n * Z[n] * 2 * binomial(n, r) / n
Bser = {k: sp.expand(v) for k, v in Bser.items() if v != 0}
Hser = bs_exp(Bser)
# h[k][i] = [eta^i] H^{(2k)}(eta/2, 0)
h = {}
for k in range(NMAX//2+1):
    for i in range(0, NMAX+1):
        v = Hser.get((i, 2*k), 0)
        if v != 0: h[(k,i)] = sp.expand(factorial(2*k)*v / 2**i)
S = {}
for k in range(NMAX//2+1):
    for i in range(0, NMAX+1):
        if (k,i) in h: S[2*k+i] = S.get(2*k+i, 0) + c[k]*h[(k,i)]
S = ts_mul(K, S)

# a_N
eta_d = {n: (-1)**n for n in range(1, NMAX+1)}
etaPow = {0: {0: sp.Integer(1)}}
for n in range(1, NMAX+1): etaPow[n] = ts_mul(etaPow[n-1], eta_d)
Sd = {0: sp.expand(S.get(0,0))}
for n in range(1, NMAX+1):
    if n in S:
        for j, v in etaPow[n].items(): Sd[j] = Sd.get(j, 0) + S[n]*v
a = {n: sp.expand(4*Sd.get(n, 0)) for n in range(NMAX+1)}

# lambda_{n,j} = coeff of Z_n^1 in s_j ; [Z_n] a_N direct
lam = {}
for n in range(3, 8, 2):
    for j in range(n, NMAX+1):
        lam[(n,j)] = sp.expand(S.get(j,0)).coeff(Z[n], 1)

print("=== Table: [Z_n] a_N (direct extraction) ===")
for n in [3, 5, 7]:
    print(f"\nn = {n}:")
    for N in range(n, 13):
        cf = sp.expand(a[N]).coeff(Z[n], 1)
        print(f"  N={N:>2}: {cf}")

print("\n=== Verify binomial transform: [Z_n]a_N = 4(-1)^N sum_j binom(N-1,j-1) lam_{n,j} ===")
ok = True
for n in [3,5,7]:
    for N in range(n, 13):
        direct = sp.expand(a[N]).coeff(Z[n], 1)
        bt = 4*(-1)**N * sum(binomial(N-1, j-1)*lam[(n,j)] for j in range(n, N+1))
        if sp.simplify(direct - bt) != 0: ok = False; print("MISMATCH", n, N)
print("binomial transform verified:", ok)

print("\n=== Verify constant part: [Z_n]a_N|_{other zetas=0} = kappa_n (-1)^{n+N} binom(N-1,n-1) ===")
kappa = {3: Rational(1), 5: Rational(-3,4), 7: Rational(129,16)}
ok2 = True
for n in [3,5,7]:
    for N in range(n, 13):
        cf = sp.expand(a[N]).coeff(Z[n], 1)
        const = sp.expand(cf.subs({Z[m]: 0 for m in range(2, NMAX+1) if m != n}))
        pred = kappa[n]*(-1)**(n+N)*binomial(N-1, n-1)
        if sp.simplify(const - pred) != 0: ok2 = False; print("MISMATCH", n, N, const, pred)
print("constant-part formula verified:", ok2)

# Verify L_n(eta) closed formula for lambda: lam_{n,j} = s_{j-n}/n - 2 sum_k c_k sum_{r even<=min(2k,n)} binom(2k,r)(n-1)!/(n-r)! 2^{r-n} h[(k-r/2, j-2k-n+r)]
def lam_formula(n, j):
    val = S.get(j-n, 0)/n
    for k in range(0, NMAX//2+1):
        for r in range(0, min(2*k, n)+1, 2):
            i = j - 2*k - n + r
            kk = k - r//2
            if i >= 0 and (kk, i) in h:
                val -= 2*c[k]*binomial(2*k, r)*Rational(factorial(n-1), factorial(n-r))*Rational(2)**(r-n)*h[(kk,i)]
    return sp.expand(val)
ok3 = all(sp.simplify(lam_formula(n,j) - lam[(n,j)])==0 for n in [3,5,7] for j in range(n, NMAX+1))
print("L_n(eta) closed formula for lambda_{n,j} verified:", ok3)
