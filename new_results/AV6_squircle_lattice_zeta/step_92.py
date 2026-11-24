code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, log, power, fsum, diff as mpdiff
setdps(60)
s = 4
def zetap(a): return mpdiff(lambda x: zeta(x), a)
# A and B by shells to K with tails
def shell_sums(K):
    A = mpf(0); B = mpf(0)
    for k in range(1, K+1):
        sig = fsum([m0*log(m0) for m0 in range(1, k)])
        A += 4*k*power(k,-s)*log(k)
        B += power(k,-s-1)*(4*sig + 2*k*log(k))
    # tails: A: 4 sum k^{1-s} ln k -> -4 zeta'(s-1) tail; B: 4 sum k^{-s-1} sigma + 2 sum k^{-s} ln k
    def HTln(sig):
        return -zetap(sig) + fsum([power(m,-sig)*log(m) for m in range(1,K+1)])
    A += 4*HTln(s-1)
    # B tail: 4 sum_{k>K} k^{-s-1} sigma(k) + 2 HTln(s)
    # sigma(k) ~ (k^2/2)ln(k-1) - (k-1)^2/4 + ... : use leading: (k^2/2 ln k - k^2/4 - (k ln k)/2 + ln k /12 + lnA)
    import mpmath
    lnA = log(mpmath.glaisher)
    def HT(sig):
        return zeta(sig) - fsum([power(m,-sig) for m in range(1,K+1)])
    Bt = mpf(0)
    Bt += 4*( mpf(1)/2*HTln(s-1) - mpf(1)/4*HT(s-1) - mpf(1)/2*HTln(s) + mpf(1)/12*HTln(s+1) + lnA*HT(s+1) )
    # correction from sigma(k) vs leading: delta(k) = sigma(k) - leading ~ O(ln k /k): sum k^{-s-1} delta:
    # delta(k) = -1/(12)ln k? no: sigma(k) = lnH(k-1): exact leading: (k^2/2-k/2+1/12)ln(k-1) - (k-1)^2/4 + lnA
    # our Bt used: (k^2/2)ln k - k^2/4 - (k ln k)/2 + (ln k)/12 + lnA : the difference sigma - this:
    # = (k^2/2-k/2+1/12)[ln(k-1)-ln k] + (k^2-k)/2·... : ln(k-1)-ln k = -sum_{i>=1} 1/(i k^i)
    for i in range(1, 40):
        Bt += 4*( -mpf(1)/(2*i)*(HT(s-1+i)) + mpf(1)/(2*i)*HT(s+i) - mpf(1)/(12*i)*HT(s+1+i) )
    Bt += 2*HTln(s)
    B += Bt
    return s*A, s*B
A_, B_ = shell_sums(2000)
print("A =", nstr(A_, 25))
print("B =", nstr(B_, 25))
print("A - B = c1 =", nstr(A_-B_, 25))
print("A expected = -4s zeta'(3) =", nstr(-4*s*zetap(s-1), 25))
print("4 zeta'(3) =", nstr(4*zetap(s-1), 25))
'''
exec(code)
