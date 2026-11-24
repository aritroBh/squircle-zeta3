from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, ci, zeta, euler, psi)
mp.dps = 100
def N(x,d=40): return nstr(x,d)
beta = lambda s: (zeta(s, mpf(1)/4) - zeta(s, mpf(3)/4))/4**s
def betap(s):  return -nsum(lambda k: (-1)**k*log(2*k+1)/(2*k+1)**s, [0, inf])
G = beta(2); Lp = log(pi)
def bp2(n):
    Jc = -( (ci(n*pi)-euler-((-1)**n)*log(pi)-log(n))/n + (1-(-1)**n)/n )/n
    Kc = ((-1)**n - 1)/n**2
    return (2/pi)*( -pi**2*log(pi)*(-1)**n/n + (2*Jc + Kc)/n )
def bp2_quad(n): return (2/pi)*quad(lambda x: x**2*log(x)*sin(n*x), [0,pi])

# check bp2 formula at several n against quadrature
for n in [1,3,5,7]:
    print(n, nstr(bp2(n),25), nstr(bp2_quad(n),25), nstr(abs(bp2(n)-bp2_quad(n)),3))

Tdir = nsum(lambda k: bp2(2*k+1)*(-1)**k, [0, inf])
print("Tdir =", N(Tdir,40))
print("exact=", N((pi/2)**2*log(pi/2),40))
