code = r'''
"""Exact u_{q,r}(s) for U_q(p) = p^{-(s-1)} sum_{r>=0} u_{q,r}(s) p^{-r}.
Machinery: u_{q,r}(s) = Gamma(s-1) zeta(s-1) sum_{l>=1} (-1)^{l+1}/l eta_l(s-1) e_{r,l}(s)  [for q=1; generalized for q],
with eta_l(x) = (-1)^{l-1}/(l-1)! * (E^{-1}-1)_{l-1} eta(x)  (falling factorial in the shift operator),
e_{r,l} = [x^r] ( exp( sum_{m>=2} q_m^{(l)} x^{m-1} ) - 1 )^l,  q_m^{(l)} = -(s-1)_m/m * eta_l(s+m-1)/eta_l(s-1).
For general q>=1: the ln^q expansion brings q! |s(n,q)|-type coefficients; we handle q=1 first.
"""
from mpmath import mp, mpf, zeta, gamma, fac, binomial, log, power, nstr

def eta(x):
    return (1-power(2, 1-x))*zeta(x)

def stirling_signed(m, t):
    # signed Stirling first kind s(m,t)
    S = [[0]*(t+1) for _ in range(m+1)]
    S[0][0] = 1
    for a in range(1, m+1):
        for b in range(0, min(a, t)+1):
            S[a][b] = S[a-1][b-1] - (a-1)*S[a-1][b] if b > 0 else -(a-1)*S[a-1][b]
    return S[m][t]

def eta_l(l, x):
    """eta_l(x) = sum_{j>=l} (-1)^{j+l} binom(j-1,l-1) j^{-x} = (-1)^{l-1}/(l-1)! (E^{-1}-1)_{l-1} eta(x)."""
    if l == 1:
        return eta(x)
    m = l-1
    # falling factorial (y)_m = sum_t s(m,t) y^t at y = E^{-1}-1
    # expand (E^{-1}-1)^t eta(x) = sum_{r=0}^t binom(t,r)(-1)^{t-r} eta(x-r)
    tot = mpf(0)
    for t in range(m+1):
        st = stirling_signed(m, t)
        if st == 0: continue
        inner = mpf(0)
        for r in range(t+1):
            inner += binomial(t, r)*((-1)**(t-r))*eta(x-r)
        tot += st*inner
    return ((-1)**m)/fac(m) * tot

def q_m(l, m, s):
    """q_m^{(l)} = -(s-1)_m/m * eta_l(s+m-1)/eta_l(s-1)."""
    poch = gamma(s-1+m)/gamma(s-1)
    return -poch/m * eta_l(l, s+m-1)/eta_l(l, s-1)

def E_series(l, s, R):
    """E^{(l)}(x) = exp(sum_{m>=2} q_m x^{m-1}) - 1 as series to x^R."""
    qs = [mpf(0)] + [q_m(l, m, s) for m in range(2, R+2)]  # qs[i] = q_{i+1}
    # exp of the series g(x) = sum_{i>=1} qs[i] x^i
    g = qs[:]  # g[0]=0
    E = [mpf(0)]*(R+1)
    E[0] = 1
    for n in range(1, R+1):
        acc = mpf(0)
        for k in range(1, n+1):
            acc += k*g[k]*E[n-k]
        E[n] = acc/n
    # subtract 1: exp(g) - 1 -> drop constant
    E[0] = 0
    return E  # E[i] = coeff of x^i

def e_rl(l, s, r):
    """[x^r] E^{(l)}(x)^l."""
    E = E_series(l, s, r)
    # power E^l
    P = [mpf(0)]*(r+1); P[0] = 1
    for _ in range(l):
        Q = [mpf(0)]*(r+1)
        for i in range(r+1):
            if P[i] == 0: continue
            for j in range(r+1-i):
                if E[j] != 0: Q[i+j] += P[i]*E[j]
        P = Q
    return P[r]

def u_1_r(s, r, lmax=None):
    """u_{1,r}(s): U_1(p) = p^{-(s-1)} sum_r u_{1,r} p^{-r}."""
    s = mpf(s)
    if lmax is None: lmax = r+1
    tot = mpf(0)
    for l in range(1, lmax+1):
        tot += ((-1)**(l+1))/l * eta_l(l, s-1) * e_rl(l, s, r)
    return gamma(s-1)*zeta(s-1)*tot
'''
with open("umach.py","w") as f: f.write(code)
test = r'''
from umach import *
from squircle import setdps
from mpmath import mp, nstr, zeta, gamma, power
setdps(80)
s = 4
# checks
print("eta_l(2, 3) =", nstr(eta_l(2, 3), 20), " should = eta(2)-eta(3) =", nstr(eta(2)-eta(3), 20))
print("eta_l(3, 5) =", nstr(eta_l(3, 5), 20))
u0 = gamma(s-1)*zeta(s-1)*eta(s)
print("u_{1,0} =", nstr(u0, 25), " (continuum)")
for r in range(1, 6):
    print(f"u_{{1,{r}}} =", nstr(u_1_r(s, r), 25))
'''
exec(test)
