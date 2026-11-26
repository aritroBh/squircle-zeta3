"""r01 (referee, own code): exact formal Taylor coefficients mu_{k,n} of M_k(alpha) = sum_n mu_{k,n} c^{2n}.

Derivation (from the paper's definition, written independently of g8_allN / g6_common):
  length of x e1 + y e2, x = t^u, y = (1-t)^u:  2 + 2u M(u),  M = 1/2 int_0^1 (sqrt(a^2+b^2-2cos(alpha) ab) - a - b) dt,
  a = t^(u-1), b = (1-t)^(u-1);  a^2+b^2-2cos(alpha)ab = (a+b)^2 - 4c^2 ab  (cos alpha = 2c^2-1).
  => mu_n(u) = 1/2 binom(1/2,n) (-4)^n J_n(u),  J_n(u) = int_0^1 (ab)^n (a+b)^(1-2n) dt.
  With g = t/(1-t):  J_n = int_0^oo g^{(n-1)(1-u)} (1+g)^(-1-u) (1+g^(1-u))^(1-2n) dg.
  (1+g^(1-u))^(1-2n) = sum_i binom(1-2n,i)(1+g)^(1-2n-i) g^i (g^-u - 1)^i ; the i-th term is O(u^i), so i <= k suffices
  for [u^k], and (g^-u-1)^i = sum_j binom(i,j)(-1)^(i-j) g^(-ju).  Each piece is a Beta integral:
  K_{i,j}(u) = Gamma(n+i+u(1-n-j)) Gamma(n+u(n+j)) / Gamma(2n+i+u).
  log Gamma(N+eps) = logGamma(N) + psi(N) eps + sum_{m>=2} psi^(m-1)(N) eps^m/m!,
  psi(N) = -gamma + H_{N-1} (gamma cancels: (1-n-j)+(n+j)-1 = 0),
  psi^(m-1)(N) = (-1)^m (m-1)! (zeta(m) - H^{(m)}_{N-1})   -- zeta(m) kept as a formal symbol.
Output: mu[k][n] = {zeta monomial (e2,e3,e4,e5,...): Fraction}."""
from fractions import Fraction as Fr
from math import comb, factorial
import pickle, sys

import os
KMAX = int(os.environ.get('KMAX', '5'))
NZ = KMAX - 1          # symbols zeta(2..KMAX)

def zadd(A, B, s=Fr(1)):
    C = dict(A)
    for m, v in B.items():
        C[m] = C.get(m, 0) + s*v
        if C[m] == 0: del C[m]
    return C

def zmul(A, B):
    C = {}
    for m1, v1 in A.items():
        for m2, v2 in B.items():
            m = tuple(x+y for x, y in zip(m1, m2))
            if sum((i+2)*e for i, e in enumerate(m)) > KMAX: continue   # weight > KMAX never needed
            C[m] = C.get(m, 0) + v1*v2
    return {m: v for m, v in C.items() if v != 0}

ONEM = (0,)*NZ
def zconst(q): return {ONEM: Fr(q)} if q else {}
def zsym(m):
    t = [0]*NZ; t[m-2] = 1; return {tuple(t): Fr(1)}

def ser_mul(A, B):
    C = [dict() for _ in range(KMAX+1)]
    for i in range(KMAX+1):
        if not A[i]: continue
        for j in range(KMAX+1-i):
            if B[j]: C[i+j] = zadd(C[i+j], zmul(A[i], B[j]))
    return C

def ser_exp(L):           # L[0] must be 0
    assert not L[0]
    E = [zconst(1)] + [dict() for _ in range(KMAX)]
    P = [zconst(1)] + [dict() for _ in range(KMAX)]
    for t in range(1, KMAX+1):
        P = ser_mul(P, L)
        P = [{m: v/t for m, v in X.items()} for X in P]
        E = [zadd(E[i], P[i]) for i in range(KMAX+1)]
    return E

_H = {}
def Hgen(N, r):           # H^{(r)}_{N} = sum_{i=1}^N 1/i^r
    key = (N, r)
    if key not in _H:
        _H[key] = sum((Fr(1, i**r) for i in range(1, N+1)), Fr(0))
    return _H[key]

def psi_coef(N, m):       # psi^{(m-1)}(N) without -gamma, as zeta-poly
    if m == 1: return zconst(Hgen(N-1, 1))
    f = (-1)**m * factorial(m-1)
    return zadd(zsym(m) if m <= KMAX else {}, zconst(-Hgen(N-1, m)))  if f == 1 else \
           {mm: f*v for mm, v in zadd(zsym(m) if m <= KMAX else {}, zconst(-Hgen(N-1, m))).items()}

def gamma_ratio_series(args):
    """args: list of (N, a, sign): prod Gamma(N + a u)^sign, normalised by prod Gamma(N)^sign; returns series of exp(L)."""
    L = [dict() for _ in range(KMAX+1)]
    for N, a, sg in args:
        for m in range(1, KMAX+1):
            pc = psi_coef(N, m)
            f = Fr(a)**m / factorial(m) * sg
            L[m] = zadd(L[m], {mm: f*v for mm, v in pc.items()})
    return ser_exp(L)

def gamma_int(N):
    return factorial(N-1)

def binom_gen(x, i):       # generalized binomial, x a Fraction
    r = Fr(1)
    for t in range(i): r = r*(x - t)/(t+1)
    return r

def mu(n, k):
    """{zeta monomial: Fraction} = [u^k] mu_n(u)."""
    tot = {}
    for i in range(0, k+1):
        bi = binom_gen(Fr(1-2*n), i)
        for j in range(0, i+1):
            cf = bi*comb(i, j)*(-1)**(i-j)
            base = Fr(gamma_int(n+i)*gamma_int(n), gamma_int(2*n+i))
            S = gamma_ratio_series([(n+i, 1-n-j, 1), (n, n+j, 1), (2*n+i, 1, -1)])
            tot = zadd(tot, {m: cf*base*v for m, v in S[k].items()})
    pre = Fr(1, 2)*binom_gen(Fr(1, 2), n)*(-4)**n
    return {m: pre*v for m, v in tot.items()}

def name(m):
    s = "*".join(f"z{i+2}" + (f"^{e}" if e > 1 else "") for i, e in enumerate(m) if e)
    return s or "1"

if __name__ == '__main__':
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    # sanity: mu_{0,n} = -1/(n(2n-1)); [z2]mu_{2,n} = 1/(4n^2-1); [z3]mu_{3,n} = -1/(4n^2-1)
    ok0 = all(mu(n, 0) == {ONEM: Fr(-1, n*(2*n-1))} for n in range(1, 30))
    z2 = zsym(2); z3 = zsym(3)
    ok2 = all(mu(n, 2).get(list(z2)[0], 0) == Fr(1, 4*n*n-1) for n in range(1, 30))
    ok3 = all(mu(n, 3).get(list(z3)[0], 0) == Fr(-1, 4*n*n-1) for n in range(1, 30))
    print("sanity mu_{0,n} = -1/(n(2n-1)), n<30:", ok0)
    print("sanity [z2]mu_{2,n} = 1/(4n^2-1), n<30 (Thm flat):", ok2)
    print("sanity [z3]mu_{3,n} = -1/(4n^2-1), n<30 (Thm flat):", ok3)
    out = {k: {} for k in range(0, KMAX+1)}
    for k in range(0, KMAX+1):
        for n in range(1, NMAX+1):
            out[k][n] = mu(n, k)
        print(f"k={k}: monomials occurring for n<={NMAX}:", sorted({name(m) for n in out[k] for m in out[k][n]}), flush=True)
    for n in (1, 2, 3):
        print(f"mu_5,{n} =", {name(m): str(v) for m, v in sorted(out[5][n].items())})
    pickle.dump(out, open(f"r01_mu_k0to5_n{NMAX}.pkl", "wb"))
