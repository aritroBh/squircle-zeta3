
import mpmath as mp
mp.mp.dps = 55

def S12(z):
    ell = mp.log(1-z)
    return mp.zeta(3) - mp.polylog(3, 1-z) + ell*mp.polylog(2, 1-z) + mp.mpf('0.5')*ell**2*mp.log(z)
def cumulants(th):
    z = mp.sin(th/2)**2; w = 1/mp.tan(th/2)**2; ell = mp.log(mp.cos(th/2)**2); u = w*ell
    return (1+u, 1+mp.pi**2/3-2*mp.polylog(2,z)-u**2, -12*w*S12(z)+2*(1+u**3))
def hc(s, th):
    c = mp.cos(th)
    return 1 + s - mp.sqrt(1 + 2*c*s + s**2)
def Jvals(th, a, A, B, n):
    k1, k2, k3 = cumulants(th)
    Xi = mp.log(n*a) + k1
    z = mp.sin(th/2)**2; c = mp.cos(th)
    mu1 = Xi; mu2 = Xi**2 + k2; mu3 = Xi**3 + 3*k2*Xi + k3
    J1 = 2*z*mu1/a
    J2 = z*(A*(mu2 - 2*mu1) - a**2*mu2 + 2*a*(1-z))/a**3
    J3 = ( z*( a**4*mu3 - (3*A*a**2 + B*a)*(mu3 - 3*mu2) + 3*A**2*(mu3 - 4*mu2 + 2*mu1) )/(3*a**5)
           - 4*A*z*(1-z)/a**4 - c*z*(1-z)/a**3 )
    return J1, J2, J3
def predict(th, a, A, B, n):
    eps = 1/(n-1); J1,J2,J3 = Jvals(th, a, A, B, mp.mpf(n))
    return J1*eps + J2*eps**2 + J3*eps**3
def quad(f, fp, th, n):
    nn = mp.mpf(n)
    def yprime(x):
        fx = f(x)
        if fx <= 0: return mp.mpf(0)
        return nn*mp.exp((nn-1)*mp.log(fx))*fp(x)
    def integrand(v):
        t = mp.exp(v)/nn
        return hc(yprime(1-t), th)*mp.exp(v)/nn
    return mp.quad(integrand, [-130, 0, 3, mp.log(nn)], maxdegree=8)

# f1: (x+x^2)/2: a=3/2, A=1, B=0
f1 = lambda x: (x+x**2)/2; f1p = lambda x: (1+2*x)/2
# f2: (3x - x^3)/2: f(0)=0,f(1)=1, a = (3-3)/2=0 NO. use f2 = x + (x-1)x^2/5 hmm
# f2: x + 0.4 x^2 (x-1): f(1)=1, increasing? f'= 1+0.8x-1.2x^2: at x=1: 0.6>0; min at... f' >= 1-1.2x^2... at x=1 =0.6 ok increasing
f2 = lambda x: x + mp.mpf('0.4')*x**2*(x-1)
f2p = lambda x: 1 + mp.mpf('0.8')*x - mp.mpf('1.2')*x**2
# a=0.6, A=f''(1)= (0.8-2.4x)|_1 = -1.6, B=f'''(1) = -2.4
cases = [('quad', f1, f1p, (mp.mpf('1.5'), mp.mpf('1'), mp.mpf('0'))),
         ('cubic', f2, f2p, (mp.mpf('0.6'), mp.mpf('-1.6'), mp.mpf('-2.4')))]
for name, f, fp, (a,A,B) in cases:
    for nm, th in [('pi/3', mp.pi/3), ('pi/2', mp.pi/2), ('2pi/3', 2*mp.pi/3)]:
        for n in [1000, 100000]:
            num = quad(f, fp, th, n); pre = predict(th, a, A, B, n)
            eps = 1/mp.mpf(n-1); sc = eps**4*mp.log(n)**4
            print(f" f={name} th={nm:6s} n={n:6d}: 2-L={mp.nstr(num,16)} resid={mp.nstr(num-pre,3)} scaled={mp.nstr((num-pre)/sc,4)}")
