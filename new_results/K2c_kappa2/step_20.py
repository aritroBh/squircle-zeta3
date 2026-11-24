
import mpmath
from mpmath import mp, mpf, mpc, sqrt, log, exp, pi, gamma
mp.dps = 70
exec(open('<project>/research/08_theory/K2c_kappa2/step_19.py').read().split("# sanity")[0])

def beta1_formula(a):
    qq = mp.sqrt(1j*a)
    return 3*g1f(a)*qq/2 - qq**3/12 + 5*qq/4 + 1/(2*qq)
def beta2_formula(a):
    qq = mp.sqrt(1j*a); g1=g1f(a); g2=g2f(a)
    return (-g1*qq**4/8 + 21*g1*qq**2/8 + 9*g1/4 + 15*g2*qq**2/4
            + qq**6/288 - 5*qq**4/48 + 65*qq**2/96 + mp.mpf(15)/8)

a = mp.pi/2; A = mp.sqrt(2*a*mp.sin(a))
b1 = beta1_formula(a); b2 = beta2_formula(a)
for N in [1500, 3200, 10000]:
    th = mp.sqrt(2*a*N)*(1+1j)
    sig = 1/mp.sqrt(mpf(N))
    P = -(2/A)*N*mp.e**th*Iplus(N, a)
    print('N=%d  P = %s' % (N, mp.nstr(P, 18)))
    print('        (P-1-b1 s)/s^2 = %s   (formula b2 = %s)' % (mp.nstr((P-1-b1*sig)/sig**2, 18), mp.nstr(b2, 8)))
# also Jm/model for m=1,2 at N=3200
def saddle_xs(N, Z):
    return (-Z + sqrt(Z**2 + 4*N*Z))/2
N=3200; Z=1j*a; x1=saddle_xs(N,Z)
d2 = N*(1/x1**2 - 1/(x1+Z)**2); phis = x1 + N*log(1+Z/x1)
for m in [0,1,2]:
    model = exp(-phis)*sqrt(2*pi/d2)*(x1+Z)**(m+mpf(3)/2)/x1
    print('m=%d Jm/model = %s' % (m, mp.nstr(Jm(N,Z,m)/model, 10)))
