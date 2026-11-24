
import mpmath as mp, time
mp.mp.dps = 55

def hc(s, th):
    c = mp.cos(th)
    return 1 + s - mp.sqrt(1 + 2*c*s + s**2)

def quad_deficit_sin(th, n):
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
    return mp.quad(integrand, [-80, mid-15, mid+5, mp.log(nn)], maxdegree=8)

# contact law: -log f(1-tau) ~ beta tau^m, m=2, beta=pi^2/8
beta = mp.pi**2/8; m = 2
def uV_solve(V):
    # u^2 + log u = log V, u>0
    return mp.findroot(lambda u: u**2 + mp.log(u) - mp.log(V), 1.2)

t0=time.time()
print("f = sin(pi x/2),  f'(1)=0 -> higher contact m=2, beta=pi^2/8")
for name, th in [('pi/3', mp.pi/3), ('pi/2', mp.pi/2), ('2pi/3', 2*mp.pi/3)]:
    c = mp.cos(th)
    for n in [100, 1000, 10000, 100000]:
        num = quad_deficit_sin(th, n)
        nn = mp.mpf(n)
        leading = (1-c)*mp.sqrt(mp.log(nn)/(m*beta*nn))
        V = 2*mp.sqrt(beta*nn)
        uV = uV_solve(V)
        refined = (1-c)*(uV + mp.mpf('0.5'))/mp.sqrt(beta*nn)
        print(f" th={name:6s} n={n:6d}: 2-L={mp.nstr(num,15)} ratio_to_leading={mp.nstr(num/leading,8)}  num-refined={mp.nstr(num-refined,4)}  rel={mp.nstr(abs(num-refined)/num,3)}")
print('elapsed', time.time()-t0)
