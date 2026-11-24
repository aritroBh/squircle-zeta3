
import mpmath
from mpmath import mp, mpf, mpc, sqrt, log, exp, pi, gamma
mp.dps = 60

def gam(a): return mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4 - a/2))
def g1f(a): return -mp.mpf(3)/4 - 1j*mp.e**(1j*a)/(4*mp.sin(a))
def g2f(a):
    s = mp.sin(a); E = mp.e**(1j*a)
    return mp.mpf(29)/96 + E**2*(mp.mpf(1)/2 + 1/(32*s**2)) + (1j*E/s)*(E**2/4 - mp.mpf(3)/16)
def saddle_xs(N, Z):
    return (-Z + sqrt(Z**2 + 4*N*Z))/2

def Jm(N, Z, m):
    x1 = saddle_xs(N, Z)
    def f(x):
        return exp(-x)*(1+Z/x)**(-N)*(x+Z)**(m+mpf(3)/2)/x
    I1 = x1*mp.quadgl(lambda t: f(t*x1), [0, mpf('0.2'), mpf('0.5'), mpf('0.8'), 1])
    d2 = N*(1/x1**2 - 1/(x1+Z)**2); ad2 = abs(d2)
    v = exp(-mpc(0,1)*mp.arg(d2)/2)
    T = 11.0/sqrt(ad2)
    I2 = (T*v)*mp.quadgl(lambda t: f(x1+t*T*v), [0, mpf('0.3'), mpf('0.7'), 1])
    x2 = x1+T*v
    I3 = mp.quadgl(lambda s: f(x2+s), [0, mpf('5'), mpf('20'), mp.inf])
    return I1+I2+I3

def Iplus(N, a):
    Z = 1j*a; tot = 0
    gms = [mp.mpf(1), g1f(a), g2f(a)]
    for m in range(3):
        bf = mp.gamma(m+mp.mpf('1.5'))*mp.gamma(N-m-mp.mpf('1.5'))/mp.gamma(N)
        tot += gms[m]*bf*Jm(N,Z,m)
    return (+1j/mp.pi)*gam(a)*tot

a = mp.pi/2; A = mp.sqrt(2*a*mp.sin(a))
b1 = 3*g1f(a)*mp.sqrt(1j*a)/2 - mp.sqrt(1j*a)**3/12 + 5*mp.sqrt(1j*a)/4 + 1/(2*mp.sqrt(1j*a))
qq = mp.sqrt(1j*a); g1=g1f(a); g2=g2f(a)
b2 = (-g1*qq**4/8 + 21*g1*qq**2/8 + 9*g1/4 + 15*g2*qq**2/4 + qq**6/288 - 5*qq**4/48 + 65*qq**2/96 + mp.mpf(15)/8)
print('b1 =', mp.nstr(b1, 20))
print('b2 =', mp.nstr(b2, 20))
for N in [1500, 3200, 10000]:
    th = mp.sqrt(2*a*N)*(1+1j)
    sig = 1/mp.sqrt(mpf(N))
    Ip = Iplus(N, a)
    P = -(2/A)*N*mp.e**th*Ip
    pred = 1 + b1*sig + b2*sig**2
    print('N=%d sig=%.5f' % (N, float(sig)))
    print('  P    =', mp.nstr(P, 20))
    print('  pred =', mp.nstr(pred, 20))
    print('  (P-1-b1*sig)/sig^2 =', mp.nstr((P-1-b1*sig)/sig**2, 20))
