
"""CORRECT u_{1,r}(s) via Bell-monomial reduction.
E_{r,l}(s) = sum_{(a): sum a_m(m-1)=r} l! S(A,l) [prod_m (-1/m)^{a_m}/a_m!] (s-1)_{r+A} eta_l(s+r+A-1)/eta_l(s-1)
u_{1,r}(s) = Gamma(s-1) zeta(s-1) sum_{l=1}^{r+1} (-1)^{l+1}/l eta_l(s-1) E_{r,l}(s)
"""
from mpmath import mp, mpf, zeta, gamma, fac, binomial, log, power, nstr

def eta(x):
    from mpmath import log as mlog
    if abs(x-1) < mpf('1e-60'): return mlog(2)
    return (1-power(2, 1-x))*zeta(x)

def stirling1s(m, t):
    S = [[0]*(t+1) for _ in range(m+1)]
    S[0][0] = 1
    for a in range(1, m+1):
        for b in range(0, min(a, t)+1):
            S[a][b] = (S[a-1][b-1] - (a-1)*S[a-1][b]) if b > 0 else -(a-1)*S[a-1][b]
    return S[m][t]

def eta_l(l, x):
    if l == 1: return eta(x)
    m = l-1
    tot = mpf(0)
    for t in range(m+1):
        st = stirling1s(m, t)
        if st == 0: continue
        inner = mpf(0)
        for r in range(t+1):
            inner += binomial(t, r)*((-1)**(t-r))*eta(x-r)
        tot += st*inner
    return ((-1)**m)/fac(m) * tot

def stirling2(A, l):
    if A < l: return 0
    S = [[0]*(l+1) for _ in range(A+1)]
    S[0][0] = 1
    for a in range(1, A+1):
        for b in range(1, min(a, l)+1):
            S[a][b] = S[a-1][b-1] + b*S[a-1][b]
    return S[A][l]

def monomials(r):
    """all (a) = (a_2, a_3, ...) with sum a_m (m-1) = r; return list of (A, w, alpha) with
       A = sum a_m, w = prod (-1/m)^{a_m}/a_m!, alpha = r+A."""
    out = []
    def rec(m, rem, cur):
        if rem == 0:
            A = sum(cur.values())
            w = mpf(1)
            for mm, am in cur.items():
                w *= ((-1)/mm)**am / fac(am)
            out.append((A, w, r+A))
            return
        if m-1 > rem: return
        a = 0
        while a*(m-1) <= rem:
            cur[m] = a
            rec(m+1, rem - a*(m-1), cur)
            a += 1
        cur.pop(m, None)
    rec(2, r, {})
    return out

def E_rl(r, l, s):
    tot = mpf(0)
    for (A, w, alpha) in monomials(r):
        S = stirling2(A, l)
        if S == 0: continue
        poch = gamma(s-1+alpha)/gamma(s-1)
        tot += fac(l)*S*w*poch*eta_l(l, s+alpha-1)/eta_l(l, s-1)
    return tot

def u_1_r(s, r):
    s = mpf(s)
    tot = mpf(0)
    for l in range(1, r+2):
        tot += ((-1)**(l+1))/l * eta_l(l, s-1) * E_rl(r, l, s)
    return gamma(s-1)*zeta(s-1)*tot
