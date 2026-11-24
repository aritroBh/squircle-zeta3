code = r'''
"""c_1(s) = -8s sum_{k>=2} k^{-s-1} sigma(k), sigma(k) = sum_{m=1}^{k-1} m ln m = ln H(k-1).
Exact evaluation: direct sum + EM tail using ln H(n) = (n^2/2+n/2+1/12)ln n - n^2/4 + ln A + sum B_{2j+2}/(2j(2j+2) n^{2j}).
"""
from mpmath import mp, mpf, zeta, gamma, fac, binomial, log, power, nstr, fsum, loggamma, bernoulli, glaisher

def sigma_tail(K, s, J=12):
    """sum_{k>K} k^{-s-1} sigma(k), using sigma(k) = ln H(k-1) expansion at k."""
    # ln H(k-1): write in terms of k: use expansion of ln H(k-1) about k directly:
    # sigma(k) = A2(k) ln(k-1) - (k-1)^2/4 + lnA + corr(k-1), A2 = (k-1)^2/2+(k-1)/2+1/12
    # Easier: sigma(k) = (k^2/2 - k/2 + 1/12) ln(k-1) - (k-1)^2/4 + lnA + corr
    # Expand ln(k-1) = ln k + ln(1-1/k) = ln k - sum_{i>=1} 1/(i k^i)
    # We build the asymptotic series of sigma(k) in powers of 1/k with ln k terms, then sum zeta/zeta'.
    tot = mpf(0)
    # do it numerically: for k = K+1 .. K2 (large), then remainder via integral approx of the smooth tail
    K2 = K*50
    def sig(k):
        return loggamma(k) + (k - mpf(1)/2)*log(k) - k  # placeholder wrong
    return tot

def sigma_exact(k):
    return fsum([m*log(m) for m in range(1, k)])

def c1_direct(s, K=2000):
    """direct sum k=2..K + EM tail from the asymptotic form of sigma(k)."""
    P = fsum([power(k, -s-1)*sigma_exact(k) for k in range(2, K+1)])
    # tail: sigma(k) = a2(k) ln(k-1) - (k-1)^2/4 + lnA + sum_{j>=1} g_j(k)
    # with a2(k) = (k-1)^2/2 + (k-1)/2 + 1/12 = k^2/2 - k/2 + 1/12
    # and g_j from lnH expansion: ln H(n) = (n^2/2+n/2+1/12)ln n - n^2/4 + lnA + sum_{j>=1} B_{2j+2}/(2j(2j+2)) n^{-2j}
    # Write sigma(k) = S(k) + R(k), S(k) = (k^2/2 - k/2 + 1/12) ln(k-1) - (k-1)^2/4 + lnA,
    # R(k) = sum_{j>=1} B_{2j+2}/(2j(2j+2)) (k-1)^{-2j} (asymptotic).
    # sum_{k>K} k^{-s-1} S(k): expand ln(k-1) = ln k - sum_{i>=1} 1/(i k^i):
    # S(k) = (k^2/2 - k/2 + 1/12)(ln k - sum_i 1/(i k^i)) - (k^2/2 - k/2 + 1/4) + lnA
    lnA = log(glaisher)
    # collect terms: (k^2/2 - k/2 + 1/12) ln k  ->  (1/2)k^2 ln k - (1/2) k ln k + (1/12) ln k
    # -(k^2/2 - k/2 + 1/12) sum_{i>=1} 1/(i k^i) = -sum_{i>=1} (1/(2i) k^{2-i} - 1/(2i) k^{1-i} + 1/(12 i) k^{-i})
    # -(k^2/2 - k/2 + 1/4)
    # Sum over k > K of k^{-s-1} * each:
    T = mpf(0)
    HT = lambda sig: zeta(sig) - fsum([power(m, -sig) for m in range(1, K+1)])
    HTp = lambda sig: -zeta(sig, 1, derivative=1) if False else None
    # zeta'(a) via mpmath: zeta(a, derivative=1)? mpmath: zeta(s, derivative=1) not supported; use diff
    from mpmath import diff as mpdiff
    zp = lambda a: mpdiff(lambda x: zeta(x), a)
    # leading log terms:
    T += mpf(1)/2*(-zp(s-1)) - mpf(1)/2*(-zp(s)) + mpf(1)/12*(-zp(s+1))
    # wait: sum k^{-s-1} k^2 ln k = sum k^{1-s} ln k = -zeta'(s-1); k^{-s-1} k ln k = -zeta'(s); k^{-s-1} ln k = -zeta'(s+1)
    # the (1/2)k^2 ln k contributes (1/2)(-zeta'(s-1))... but HT vs full zeta: use HT tail instead:
    def HTln(sig):
        # sum_{k>K} k^{-sig} ln k
        return -mpdiff(lambda x: zeta(x), sig) + fsum([power(m, -sig)*log(m) for m in range(1, K+1)])
    T = mpf(0)
    T += mpf(1)/2*HTln(s-1) - mpf(1)/2*HTln(s) + mpf(1)/12*HTln(s+1)
    # -(k^2/2 - k/2 + 1/4): sum k^{-s-1}(-k^2/2 + k/2 - 1/4) = -(1/2)HT(s-1) + (1/2)HT(s) - (1/4)HT(s+1)
    T += -mpf(1)/2*HT(s-1) + mpf(1)/2*HT(s) - mpf(1)/4*HT(s+1)
    # + lnA * sum k^{-s-1} = lnA * HT(s+1)
    T += lnA*HT(s+1)
    # -(k^2/2 - k/2 + 1/12) sum_{i>=1} 1/(i k^i): truncate i to 40
    for i in range(1, 41):
        T += -(mpf(1)/(2*i))*HT(s-1+i) + (mpf(1)/(2*i))*HT(s+i) - (mpf(1)/(12*i))*HT(s+1+i)
    # R(k) asymptotic: sum_{j>=1} B_{2j+2}/(2j(2j+2)) (k-1)^{-2j}: expand (k-1)^{-2j} = k^{-2j}(1+2j/k + binom(2j+1,2)/k^2 + ...)
    for j in range(1, 13):
        cj = bernoulli(2*j+2)/(2*j*(2*j+2))
        # (k-1)^{-2j} = sum_{t>=0} binom(2j+t-1, t) k^{-2j-t}
        for t in range(0, 40):
            T += cj*binomial(2*j+t-1, t)*HT(s+2*j+t)
    return -8*s*(P + T)
'''
with open("c1.py","w") as f: f.write(code)
test = r'''
from c1 import *
from squircle import setdps
from mpmath import mp, nstr, zeta
setdps(70)
for s in [3,4,5,6]:
    print(f"c_1({s}) =", nstr(c1_direct(s, K=800), 45))
'''
exec(test)
with open("progress.log","a") as f: f.write("step_47: c_1(s) exact lattice sum\n")
