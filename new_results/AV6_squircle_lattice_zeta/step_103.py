code = r'''
from squircle import setdps
from mpmath import mp, nstr, log, power, fsum, mpf, zeta, diff as mpdiff
setdps(60)
s = 4; K = 1200
P = mpf(0); acc = mpf(0)
for k in range(1, K+1):
    if k >= 2: P += power(k, -s-1)*log(k)*acc
    acc += k*log(k)
print("P (fixed loop):", nstr(P, 20))
print("hand-estimate: 0.018649 + small tail")
# tail parts
def zetap(a): return mpdiff(lambda x: zeta(x), a)
def zetapp(a): return mpdiff(lambda x: zeta(x), a, 2)
def HT(sig):
    return zeta(sig) - fsum([power(m,-sig) for m in range(1, K+1)])
def HTln(sig):
    return -zetap(sig) - fsum([power(m,-sig)*log(m) for m in range(1, K+1)])
def HTln2(sig):
    return zetapp(sig) - fsum([power(m,-sig)*log(m)**2 for m in range(1, K+1)])
t1 = HTln(s+1)*acc
print("t1 = HTln(5)*acc:", nstr(t1, 15))
t2 = (1/s)*HTln2(s-1) + (1/s**2)*HTln(s-1) + mpf(1)/2*HTln2(s)
print("t2 (Hurwitz leading):", nstr(t2, 15))
'''
exec(code)
