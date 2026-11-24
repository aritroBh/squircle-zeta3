code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, quad, exp, log, power, gamma, fsum
setdps(70)
p = mpf(50); s = mpf(4); a = s/p
ap = gamma(1+1/p)
delta = mpf('2.46e-150')
term1 = 4*ap*ap*delta**(a-2/p)/(a-2/p)
term2 = delta**a/a
print("ap:", nstr(ap, 20))
print("delta**(a-2/p):", nstr(delta**(a-2/p), 20))
print("term1:", nstr(term1, 25))
print("term2:", nstr(term2, 25))
print("small:", nstr(term1-term2, 25))
# hand-check delta**0.04
print("ln delta:", nstr(log(delta), 20))
print("check e^-13.7795:", nstr(exp(mpf('-13.7795')), 20))
# ZA T1 brute check for a few m
def h(m,n): return power(power(m,p)+power(n,p),-a)
K=250; R=16
bern = [bernoulli(2*r)/fac(2*r) for r in range(R+1)]
def T1(m):
    inner = fsum([h(m, n) for n in range(1, K+1)])
    ck = taylor_h(m, K, p, a, 2*R+1)
    tail = power(m, 1-s)*Jint(mpf(K)/m, p, a) - h(m, K)/2
    for r in range(1, R+1):
        tail -= bern[r]*ck[2*r-1]*fac(2*r-1)
    return inner + tail
for m in [41, 100, 300]:
    v1 = T1(m)
    # brute: sum to 2000 with integral tail
    vb = fsum([h(m,n) for n in range(1,2001)]) + power(m,1-s)*Jint(mpf(2000)/m,p,a)
    print(f"T1({m}): {nstr(v1,30)}  brute: {nstr(vb,30)}  diff: {nstr(v1-vb,4)}")
'''
exec(code)
