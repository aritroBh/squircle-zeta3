code = r'''
from squircle import *
from c2 import G_22
from mpmath import mp, zeta, mpf, nstr, log, power, fsum, diff as mpdiff
setdps(60)
def zetapp(a): return mpdiff(lambda x: zeta(x), a, 2)
# brute G_22(4) partial
P = mpf(0); acc = mpf(0)
for k in range(1, 1501):
    if k >= 2: P += power(k, -5)*acc
    acc += k**2*log(k)**2
print("G_22(4) brute partial:", nstr(P, 22))
print("G_22(4) from code    :", nstr(G_22(4), 22))
# U_2 check at s=3: quadrant sum 8(G_22+G_x)(4) via direct sigma_22/sigma_x tails
# brute U_2 quadrant partial to K=600
U2q = mpf(0)
for k in range(1, 601):
    for m0 in range(1, k):
        a_ = m0*log(m0); b_ = (k-m0)*log(k-m0)
        U2q += 4*power(k,-5)*(a_+b_)**2
print("U_2 quadrant partial K=600:", nstr(U2q, 20))
# tail: shell = 8 sigma_22 + 8 sigma_x + 2k^2 ln^2 k; sigma_22 ~ k^3/3 ln^2 k - 2k^3/9 ln k + 2k^3/27
# sigma_x ~ k^3[ln^2 k/6 - 5 ln k/18 + 1/18]
# use these leads for tail with EM-shifts: tail_lead = sum_{k>600} k^{-5}[8(k^3/3 ln^2 k - 2k^3/9 ln k + 2k^3/27) + 8(k^3/6 ln^2k - 5k^3/18 ln k + k^3/18)]
def HT(sig):
    return zeta(sig) - fsum([power(m,-sig) for m in range(1, 601)])
def HTln(sig):
    return -mpdiff(lambda x: zeta(x), sig) - fsum([power(m,-sig)*log(m) for m in range(1, 601)])
def HTln2(sig):
    return zetapp(sig) - fsum([power(m,-sig)*log(m)**2 for m in range(1, 601)])
tail = 8*( mpf(1)/3*HTln2(2) - mpf(2)/9*HTln(2) + mpf(2)/27*HT(2)
          + mpf(1)/6*HTln2(2) - mpf(5)/18*HTln(2) + mpf(1)/18*HT(2) )
print("tail (lead approx):", nstr(tail, 15))
print("U_2 quadrant+leadtail:", nstr(U2q + tail, 20))
print("8(G_22+G_x)(4) code   :", nstr(8*(G_22(4)+__import__('c2').G_x(4)), 20))
'''
exec(code)
