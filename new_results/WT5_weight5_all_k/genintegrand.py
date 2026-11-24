
"""Closed-form M_k integrand in the Euler variable v, pure rational arithmetic.

Phi(eps) = A e^{eps l1} + B e^{eps l2} - Q' sqrt(1+u),
  u = alpha X + beta Y,  X = e^{2 eps l1}-1, Y = e^{2 eps l2}-1,
  alpha = (1-t)^2/q2, beta = t^2/q2, Q' = q/(t(1-t)), q2 = t^2+(1-t)^2, q = sqrt(q2).

[eps^k] Phi = A l1^k/k! + B l2^k/k! - Q' sum_{j=1}^k binom(1/2,j) [eps^k] u^j,
[eps^k] X^s Y^{j-s} = sum_m w(s,m) w(j-s,k-m) l1^m l2^{k-m},  w(s,m) = s! S(m,s) 2^m/m!.

So the coefficient of l1^m l2^{k-m} in the M_k integrand (-1/2 * [eps^k] Phi) is
  sqrt part: -q/(2 t(1-t) q2^k) * S_m(t),
  S_m(t) = sum_{j,s} binom(1/2,j) binom(j,s) w(s,m) w(j-s,k-m) t^{2(j-s)} (1-t)^{2s} q2^{k-j}
  plus A- term (m=k): -1/(2 k! t),   B- term (m=0): -1/(2 k! (1-t)).

In v: t=(v^2+2v-1)/(4v), 1-t=(1+2v-v^2)/(4v), q=(v^2+1)/(2 sqrt2 v), dt/dv=(v^2+1)/(4 v^2).
Coefficient * dt/dv, sqrt part:
  -sqrt2 8^k 4^{-D} v^{2k-1-D} Stilde(v) (v^2+1)^{1-2k} / ((v^2-a^2)(v^2-b^2))
  Stilde(v) = sum_s S_s (v^2+2v-1)^s (4v)^{D-s},  D = deg S.
A-part * dt/dv: -(v^2+1)/(2 k! v (v-b)(v+a));  B-part: -(v^2+1)/(2 k! v (a-v)(v+b)).
Folding i<j at the t-level: c -> c(1-t), i.e. S(t) -> S(1-t) (q2, t(1-t) invariant), A <-> B.
"""
from fractions import Fraction as Fr
from math import comb, factorial

def stirling2(n, k):
    if n < 0 or k < 0: return 0
    if n == k == 0: return 1
    if k == 0 or n == 0: return 0
    # via recurrence row-wise for small n this is fine
    prev = [0]*(k+1); prev[0] = 1
    for i in range(1, n+1):
        cur = [0]*(k+1)
        for j in range(1, min(i,k)+1):
            cur[j] = j*prev[j] + prev[j-1]
        prev = cur
    return prev[k] if k <= n else 0

_S = {}
def S2_mem(n, k):
    key = (n,k)
    if key not in _S:
        _S[key] = stirling2(n,k)
    return _S[key]

def w(s, m):
    """s! S(m,s) 2^m / m!   (convention: s=0 gives delta_{m,0})"""
    if s == 0:
        return Fr(1) if m == 0 else Fr(0)
    if m < s: return Fr(0)
    return Fr(factorial(s) * S2_mem(m, s) * 2**m, factorial(m))

_HALF = Fr(1,2)
def binom_half(j):
    r = Fr(1)
    for t in range(j):
        r *= (_HALF - t)
    return r / factorial(j)

# polynomial arithmetic over Q (lists low-to-high of Fr)
def padd(P, Q):
    n = max(len(P), len(Q)); out = [Fr(0)]*n
    for i,c in enumerate(P): out[i] += c
    for i,c in enumerate(Q): out[i] += c
    return out
def pscale(P, q):
    return [c*q for c in P]
def pmul(P, Q):
    if not P or not Q: return []
    out = [Fr(0)]*(len(P)+len(Q)-1)
    for a,ca in enumerate(P):
        if ca == 0: continue
        for b,cb in enumerate(Q):
            if cb: out[a+b] += ca*cb
    return out
def pshift(P, n):  # multiply by t^n
    return [Fr(0)]*n + P
def ptrim(P):
    while P and P[-1] == 0: P.pop()
    return P
def p_comp_1mt(P):
    # P(1-t)
    out = [Fr(0)]
    for s, c in enumerate(P):
        if c == 0: continue
        term = [Fr(1)]
        for _ in range(s):
            term = pmul(term, [Fr(1), Fr(-1)])
        out = padd(out, pscale(term, c))
    return ptrim(out)

Q2 = [Fr(1), Fr(-2), Fr(2)]    # q2 = t^2 + (1-t)^2 = 1 - 2t + 2t^2

def S_poly(k, m):
    """S_m(t) as above."""
    out = []
    for j in range(1, k+1):
        bh = binom_half(j)
        for s in range(0, j+1):
            w1 = w(s, m); w2 = w(j-s, k-m)
            if w1 == 0 or w2 == 0: continue
            c = bh * comb(j, s) * w1 * w2
            # t^{2(j-s)} (1-t)^{2s} q2^{k-j}
            term = [Fr(1)]
            term = pshift(term, 2*(j-s))
            onemt = [Fr(1)]
            for _ in range(2*s):
                onemt = pmul(onemt, [Fr(1), Fr(-1)])
            term = pmul(term, onemt)
            for _ in range(k-j):
                term = pmul(term, Q2)
            out = padd(out, pscale(term, c))
    return ptrim(out)

def coeffs_Q(k):
    """returns dict (m, k-m) -> (kind, data) for the unfolded M_k coefficient of l1^m l2^{k-m}.
       kind 'sqrt': data = S(t) (coeff = -q S/(2 t(1-t) q2^k))
       additionally m=k gets A-part, m=0 gets B-part (flags)."""
    out = {}
    for m in range(0, k+1):
        out[(m, k-m)] = S_poly(k, m)
    return out
