code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, fsum, power, inf, quad
setdps(50)
p = mpf(2); s = mpf(4); a = s/p

# Jint check vs numerical quad
for x in [mpf('0.4'), mpf('2.5')]:
    j1 = Jint(x, p, a)
    j2 = quad(lambda u: (1+u**p)**(-a), [x, inf])
    print("Jint", nstr(x,3), nstr(j1,15), nstr(j2,15), nstr(fabs(j1-j2),3))

# T1(m) check: sum_{n>=1} (m^p+n^p)^{-a}
m = 31; K = 120; R = 12
bern = [bernoulli(2*r)/fac(2*r) for r in range(R+1)]
def h(m,n): return (m**p + n**p)**(-a)
inner = fsum([h(m,n) for n in range(1,K+1)])
ck = taylor_h(m, K, p, a, 2*R+1)
tail = m*Jint(mpf(K)/m, p, a) - h(m,K)/2
for r in range(1,R+1):
    tail -= bern[r]*ck[2*r-1]*fac(2*r-1)
T1 = inner + tail
brute = fsum([h(m,n) for n in range(1,20000)]) 
print("T1   :", nstr(T1, 25))
print("brute:", nstr(brute, 25))
print("diff :", nstr(T1-brute, 5), " (brute tail ~", nstr(h(m,20000)/((s-1)*20000**(s-1))*20000**s,5),")")
# brute tail estimate: sum_{n>20000} h(m,n) ~ int_{20000} = m*Jint(20000/m,...)
bt = m*Jint(mpf(20000)/m, p, a)
print("brute+tail:", nstr(brute+bt, 25))
'''
exec(code)
