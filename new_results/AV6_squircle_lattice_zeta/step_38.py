code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, quad, exp, log, power, gamma, gammainc, fsum
import time
setdps(70)
p = mpf(50); s = mpf(4); a = s/p
ap = gamma(1+1/p); ga = gamma(a)
delta = mpf('2.46e-150'); ld = float(log(delta))
cache = {}
def mpw(m):
    v = cache.get(m)
    if v is None:
        v = power(m, p); cache[m] = v
    return v
def theta(t):
    f = mpf(0); cutoff = mpf(10)**(-85); m = 0
    while True:
        m += 1
        e = exp(-t*mpw(m)); f += e
        if e < cutoff: break
    return 1+2*f
# termwise int_delta^1 t^{a-1} 4f(t) dt = 4 sum_{m>=1} m^{-s} [gam(a, m^p) - gam(a, delta m^p)]
t0=time.time()
term4f = mpf(0)
for m in range(1, 3000):
    mp_ = mpw(m)
    g1 = gammainc(a, 0, mp_)          # = Gamma(a) - Gamma(a, m^p)
    g2 = gammainc(a, 0, delta*mp_)
    term4f += power(m, -s)*(g1-g2)
term4f *= 4
# m-tail beyond 3000: 4 sum_{m>3000} m^-s [gam(a,m^p)-gam(a,dm^p)] ~ 4 zeta-tail Gamma(a) (both gams -> Gamma(a) and 0)
# g1 -> Gamma(a), g2 -> gammainc(a,0,~0)=0  for m>3000 (m^50 huge, delta*m^50 still huge? delta*3000^50 = 2.46e-150*e^{50*8.006}=e^{-344.5+400.3}=e^{55.8} -> gammainc(a,0,big)=Gamma(a); so g1-g2->0!)
# indeed for m>~1400, delta*m^50 > 195 so g2 = Gamma(a) too -> g1-g2 ~ Gamma(a, delta m^p)-Gamma(a,m^p) ~ 0
print("term4f =", nstr(term4f, 45), time.time()-t0)
# quad of 4f alone on [delta,1] in x-space
def fonly(t):
    f = mpf(0); cutoff = mpf(10)**(-85); m = 0
    while True:
        m += 1
        e = exp(-t*mpw(m)); f += e
        if e < cutoff: break
    return f
pts = [ld]
while pts[-1] < -25: pts.append(pts[-1]+25)
pts.append(0.0)
I1f = mpf(0)
for i in range(len(pts)-1):
    I1f += quad(lambda x: exp(a*x)*4*fonly(exp(x)), [pts[i], pts[i+1]], method='gauss-legendre')
print("I1f quad =", nstr(I1f, 45), time.time()-t0)
print("diff:", nstr(I1f-term4f, 6))
'''
exec(code)
