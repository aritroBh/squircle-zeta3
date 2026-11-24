
import mpmath as mp, time
mp.mp.dps = 55

def hc(s, th):
    c = mp.cos(th)
    return 1 + s - mp.sqrt(1 + 2*c*s + s**2)

# Method 1: v-coordinate
def quad_sin_v(th, n):
    nn = mp.mpf(n)
    f = lambda x: mp.sin(mp.pi*x/2)
    fp = lambda x: mp.pi/2*mp.cos(mp.pi*x/2)
    def yprime(x):
        fx = f(x)
        if fx <= 0: return mp.mpf(0)
        return nn*mp.exp((nn-1)*mp.log(fx))*fp(x)
    def integrand(v):
        t = mp.exp(v)/nn
        return hc(yprime(1-t), th)*mp.exp(v)/nn
    mid = 0.5*mp.log(nn)
    knots = [-80, 0, mid, mp.log(nn)]
    return mp.quad(integrand, knots, maxdegree=8)

# Method 2: t = sqrt(n)(1-x), t in [0, sqrt n]
def quad_sin_t(th, n):
    nn = mp.mpf(n); sq = mp.sqrt(nn)
    def integrand(t):
        x = 1 - t/sq
        fx = mp.sin(mp.pi*x/2)
        yp = nn*mp.exp((nn-1)*mp.log(fx))*(mp.pi/2)*mp.cos(mp.pi*x/2)
        return hc(yp, th)/sq
    return mp.quad(integrand, [0, 3, 15, sq], maxdegree=8)

beta = mp.pi**2/8; m = 2
def uV_solve(V):
    return mp.findroot(lambda u: u**2 + mp.log(u) - mp.log(V), 1.2)

t0=time.time()
print("f=sin: two quadrature methods, then contact-law comparison")
for name, th in [('pi/3', mp.pi/3), ('pi/2', mp.pi/2), ('2pi/3', 2*mp.pi/3)]:
    c = mp.cos(th)
    for n in [100, 1000, 10000, 100000]:
        q1 = quad_sin_v(th, n); q2 = quad_sin_t(th, n)
        nn = mp.mpf(n)
        leading = (1-c)*mp.sqrt(mp.log(nn)/(m*beta*nn))
        V = 2*mp.sqrt(beta*nn); uV = uV_solve(V)
        refined = (1-c)*(uV + mp.mpf('0.5'))/mp.sqrt(beta*nn)
        print(f" th={name:6s} n={n:6d}: v={mp.nstr(q1,15)} t={mp.nstr(q2,15)} d={mp.nstr(abs(q1-q2),3)} ratio_lead={mp.nstr(q2/leading,7)} rel_ref={mp.nstr(abs(q2-refined)/q2,3)}")
print('elapsed', time.time()-t0)
