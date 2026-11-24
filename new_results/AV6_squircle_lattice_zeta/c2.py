
"""Glaisher-type sums and c_2(s).
G_2(s) = sum_k k^{-s-1} sigma_2(k), sigma_2(k) = sum_{m<k} m ln^2 m
G_l(s) = sum_k k^{-s-1} ln k sigma(k)
G_x(s) = sum_k k^{-s-1} sigma_x(k), sigma_x(k) = sum_{m<k} m ln m (k-m) ln(k-m)
c_2(s) = 2 s^2 zeta''(s-1) - s(3s+1) zeta''(s) + 4s(s+1)[G_2(s+1)+G_x(s+1)] - 8 s^2 G_l(s) - 4 s G_2(s)
"""
from mpmath import mp, mpf, zeta, gamma, fac, binomial, log, power, nstr, fsum, bernoulli, diff as mpdiff

def zetap(a): return mpdiff(lambda x: zeta(x), a)
def zetapp(a): return mpdiff(lambda x: zeta(x), a, 2)

def HT(K, sig):
    return zeta(sig) - fsum([power(m,-sig) for m in range(1, K+1)])
def HTln(K, sig):
    return -zetap(sig) - fsum([power(m,-sig)*log(m) for m in range(1, K+1)])
def HTln2(K, sig):
    return zetapp(sig) - fsum([power(m,-sig)*log(m)**2 for m in range(1, K+1)])

def sum_m_ln2(K):
    return fsum([m*log(m)**2 for m in range(1, K+1)])
def sum_m_ln(K):
    return fsum([m*log(m) for m in range(1, K+1)])

def G_2(s, K=2000, R=10):
    P = mpf(0)
    acc = mpf(0)
    for k in range(1, K+1):
        if k >= 2: P += power(k, -s-1)*acc
        acc += k*log(k)**2
    tail = HT(K, s+1)*sum_m_ln2(K)
    tail += (1/s)*HTln2(K, s-1) - mpf(1)/2*HTln2(K, s)
    for r in range(1, R+1):
        poch = gamma(s+1+2*r-1)/gamma(s+1)
        tail -= bernoulli(2*r)/fac(2*r)*poch*HTln2(K, s+2*r)
    return P + tail

def G_l(s, K=2000, R=10):
    # sum_k k^{-s-1} ln k sigma(k): direct + tail via m ln m * [-zeta'(s+1, m+1)]
    P = mpf(0)
    acc = mpf(0)
    for k in range(1, K+1):
        if k >= 2: P += power(k, -s-1)*log(k)*acc
        acc += k*log(k)
    # tail: sum_{m<=K} m ln m * HTln(K, s+1) + sum_{m>K} m ln m * [-zeta'(s+1, m+1)]
    # -zeta'(sig, m+1) = m^{1-sig} ln m/(sig-1) + m^{1-sig}/(sig-1)^2 + m^{-sig} ln m /2
    #   + sum_r B_{2r}/(2r)! [ (sig)_{2r-1} ln m + (sig)_{2r-1} (H_{sig+2r-1}-H_sig) ] m^{-sig-2r+1}
    tail = HTln(K, s+1)*sum_m_ln(K)
    tail += (1/s)*HTln2(K, s-1) + (1/s**2)*HTln(K, s-1) + mpf(1)/2*HTln2(K, s)
    def pochH(sig, r):
        # (sig)_{2r-1} and its derivative wrt sig: (sig)_{2r-1}*(H_{sig+2r-1}-H_sig)
        from mpmath import digamma
        P_ = gamma(sig+2*r-1)/gamma(sig)
        return P_, P_*(digamma(sig+2*r-1)-digamma(sig))
    for r in range(1, R+1):
        P_, dP = pochH(s+1, r)
        tail -= bernoulli(2*r)/fac(2*r)*(P_*HTln2(K, s+2*r) + dP*HTln(K, s+2*r))
    return P + tail

def G_x(s, K=1000):
    # sigma_x(k): direct to K2=4K, tail via fit of residual as (lnk,1/k) series from exact values
    K2 = 4*K
    P = mpf(0)
    for k in range(2, K2+1):
        sx = fsum([m*log(m)*(k-m)*log(k-m) for m in range(1, k)])
        P += power(k, -s-1)*sx
    # far tail k>K2: sigma_x(k) ~ k^3[ln^2 k/6 - 5 ln k/18 + 1/18] + O(k^2 ln^2 k)
    # fit residual R(k) = sigma_x(k) - k^3[...] from exact values at k = K2+1..K2+60
    res = []
    for k in range(K2+1, K2+61):
        sx = fsum([m*log(m)*(k-m)*log(k-m) for m in range(1, k)])
        lead = k**3*(log(k)**2/6 - 5*log(k)/18 + mpf(1)/18)
        res.append((k, sx - lead))
    # model: R(k) = k^2 ln^2 k * (c0 + c1/k + c2/k^2) 
    import mpmath
    A = mpmath.matrix(60, 3); b = mpmath.matrix(60, 1)
    for i, (k, rv) in enumerate(res):
        A[i,0] = k**2*log(k)**2; A[i,1] = k*log(k)**2; A[i,2] = log(k)**2
        b[i,0] = rv
    c = mpmath.lu_solve(A.T*A, A.T*b)
    tail = mpf(0)
    for k in range(K2+1, 20*K2+1):
        model = c[0,0]*k**2*log(k)**2 + c[1,0]*k*log(k)**2 + c[2,0]*log(k)**2
        tail += power(k, -s-1)*model
    # beyond 20K2 negligible (model ~ ln^2 k * k^{1-s})
    return P + tail

"""G_2^{(2)}(s) = sum_k k^{-s-1} sum_{m<k} m^2 ln^2 m."""
from mpmath import mp, mpf, zeta, gamma, fac, binomial, log, power, nstr, fsum, bernoulli, diff as mpdiff

def zetapp(a): return mpdiff(lambda x: zeta(x), a, 2)

def HT(K, sig):
    return zeta(sig) - fsum([power(m,-sig) for m in range(1, K+1)])
def HTln2(K, sig):
    return zetapp(sig) - fsum([power(m,-sig)*log(m)**2 for m in range(1, K+1)])

def G_22(s, K=2000, R=10):
    P = mpf(0)
    acc = mpf(0)
    for k in range(1, K+1):
        if k >= 2: P += power(k, -s-1)*acc
        acc += k**2*log(k)**2
    tail = HT(K, s+1)*acc
    tail += (1/s)*HTln2(K, s-2) - mpf(1)/2*HTln2(K, s-1)
    for r in range(1, R+1):
        poch = gamma(s+1+2*r-1)/gamma(s+1)
        tail -= bernoulli(2*r)/fac(2*r)*poch*HTln2(K, s+2*r-2)
    return P + tail
