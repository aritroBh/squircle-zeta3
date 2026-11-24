
import mpmath as mp
mp.mp.dps = 50

def hc(s, c):
    return 1 + s - mp.sqrt(1 + 2*c*s + s**2)
def rho(s, c):
    return mp.sqrt(1 + 2*c*s + s**2) - s - c
def I(V, c):
    return mp.quad(lambda u: hc(V*u*mp.exp(-u**2), c), [0, 1, 3, mp.inf], maxdegree=9)
def ustar(V):  # s(u*)=1 exactly: V u* e^{-u*^2} = 1
    return mp.findroot(lambda u: u**2 - mp.log(u) - mp.log(V), (mp.log(V))**mp.mpf('0.5'))
def Klayer(c):
    i1 = mp.quad(lambda s: hc(s, c)/s, [0, 1])
    i2 = mp.quad(lambda s: rho(s, c)/s, [1, mp.inf])
    return (i1 - i2)/2

for c in [mp.mpf(0), mp.mpf('0.5'), mp.mpf('-0.5')]:
    Ktot = Klayer(c) - (1-c**2)/4
    print(f"c={c}: predicted coeff K(c)={mp.nstr(Ktot,12)}")
    for lV in [3,4,5,6,8]:
        V = mp.mpf(10)**lV
        us = ustar(V)
        lhs = (I(V,c) - (1-c)*us)*us
        print(f"   V=1e{lV}: u*={mp.nstr(us,8)}  (I-(1-c)u*)u*={mp.nstr(lhs,10)}  pred={mp.nstr(Ktot,10)}  diff={mp.nstr(lhs-Ktot,5)}")
