code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, fsum, power
setdps(50)
p = mpf(2); s = mpf(4); a = s/p
def h(m,n): return (m**p + n**p)**(-a)
R=12; K=120
bern = [bernoulli(2*r)/fac(2*r) for r in range(R+1)]
for m in [31, 50, 100]:
    inner = fsum([h(m,n) for n in range(1,K+1)])
    ck = taylor_h(m, K, p, a, 2*R+1)
    tail = power(m,1-s)*Jint(mpf(K)/m, p, a) - h(m,K)/2
    for r in range(1,R+1):
        tail -= bern[r]*ck[2*r-1]*fac(2*r-1)
    T1 = inner+tail
    brute = fsum([h(m,n) for n in range(1,5000)]) + power(m,1-s)*Jint(mpf(5000)/m,p,a)
    print(m, nstr(T1,20), nstr(brute,20), nstr(T1-brute,4))
# also check Ap vs numerical
print("Aint:", nstr(Aint(p,a),20), nstr(quad(lambda u:(1+u**p)**(-a),[0,inf]),20))
# whole-line: sum_{n in Z} h(m,n) vs 2 A m^{1-s}
m=40
wl = h(m,0) + 2*fsum([h(m,n) for n in range(1,3000)]) + 2*power(m,1-s)*Jint(mpf(3000)/m,p,a)
print("whole line:", nstr(wl,20), nstr(2*Aint(p,a)*power(m,1-s),20), nstr(wl-2*Aint(p,a)*power(m,1-s),4))
'''
exec(code)
