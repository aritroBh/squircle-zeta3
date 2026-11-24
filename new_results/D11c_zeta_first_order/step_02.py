# step_02.py -- D11c(1): S(eta), a_n, kappa_n extraction + closed-form verification (sympy exact)
import sympy as sp
from sympy import symbols, Rational, factorial, binomial, sqrt, exp, series

NMAX = 20
Z = {n: symbols(f'Z{n}') for n in range(2, NMAX+1)}
eta, mu, s, d, t = symbols('eta mu s d t')

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

# c_k : sqrt(cosh z) = sum c_k z^{2k}
coshm1 = {2*k: Rational(1, factorial(2*k)) for k in range(1, NMAX//2+1)}
sqcosh = ts_pow(coshm1, Rational(1,2))
c = {k: sqcosh.get(2*k, 0) for k in range(NMAX//2+1)}
# K(eta) = exp(-sum (-1)^n Z_n eta^n/n)
mA = {n: (-1)**(n+1)*Z[n]/n for n in range(2, NMAX+1)}
K = ts_exp(mA)
# H(mu,s)
Bser = {}
for n in range(2, NMAX+1):
    for r in range(0, n+1, 2):
        Bser[(n-r,r)] = Bser.get((n-r,r), 0) + (-1)**n*Z[n]*2*binomial(n,r)/n
Bser = {k: sp.expand(v) for k, v in Bser.items() if v != 0}
Hser = bs_exp(Bser)
# S(eta) = K * sum_k c_k eta^{2k} H^{(2k)}(eta/2,0)  = P/(4 sqrt2)
S = {}
for k in range(NMAX//2+1):
    for i in range(0, NMAX+1):
        v = Hser.get((i,2*k), 0)
        if v != 0: S[2*k+i] = S.get(2*k+i,0) + c[k]*factorial(2*k)*v/2**i
S = ts_mul(K, S)
# a_n from eta = -d/(1+d)
eta_d = {n: (-1)**n for n in range(1, NMAX+1)}
etaPow = {0: {0: sp.Integer(1)}}
for n in range(1, NMAX+1): etaPow[n] = ts_mul(etaPow[n-1], eta_d)
Sd = {0: sp.expand(S.get(0,0))}
for n in range(1, NMAX+1):
    if n in S:
        for j, v in etaPow[n].items(): Sd[j] = Sd.get(j, 0) + S[n]*v
a = {n: sp.expand(4*Sd.get(n, 0)) for n in range(NMAX+1)}

print("verified coefficients (Z2=pi^2/6, Z4=pi^4/90, ...):")
for n in range(2, 8): print(f"  a_{n} =", a[n])

# kappa extraction + both closed forms
Rt = series(sqrt((1+exp(4*t))/2), t, 0, NMAX+1).removeO()
print(f"\n{'n':>3} {'kappa extracted':>26} {'V-form':>26} {'c_k-sum':>26} ok")
for n in range(3, 20, 2):
    kex = sp.expand(S.get(n,0)).coeff(Z[n],1)*4*(-1)**n
    Vn = factorial(n)*Rt.coeff(t,n)
    kV = sp.Rational(4*(-1)**n, n)*(1 - sp.Rational(2)**(1-n)*Vn)
    kc = sp.Rational(4*(-1)**n, n)*(1 - 2*sum(c[k]*factorial(n)/factorial(n-2*k)*sp.Rational(2)**(2*k-n)
                                              for k in range(0,(n-1)//2+1)))
    ok = sp.simplify(kex-kV)==0 and sp.simplify(kex-kc)==0
    print(f"{n:>3} {str(kex):>26} {str(sp.simplify(kV)):>26} {str(sp.simplify(kc)):>26} {ok}")
