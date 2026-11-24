
import mpmath as mp
mp.mp.dps = 50

def hc(s, c):
    return 1 + s - mp.sqrt(1 + 2*c*s + s**2)
def rho(s, c):
    return mp.sqrt(1 + 2*c*s + s**2) - s - c

def G(c):
    i1 = mp.quad(lambda s: hc(s, c)/s, [0, 1])
    i2 = mp.quad(lambda s: rho(s, c)/s, [1, mp.inf])
    return (i1 - i2)/2, i1, i2

def I(V, c):
    return mp.quad(lambda u: hc(V*u*mp.exp(-u**2), c), [0, 1, 3, mp.inf], maxdegree=9)
def uV_solve(V):
    return mp.findroot(lambda u: u**2 + mp.log(u) - mp.log(V), (mp.log(V))**mp.mpf('0.5'))

for c in [mp.mpf(0), mp.mpf('0.5'), mp.mpf('-0.5')]:
    Gc, i1, i2 = G(c)
    print(f"c={c}: G(c)={mp.nstr(Gc,14)} (i1={mp.nstr(i1,12)}, i2={mp.nstr(i2,12)})  G/(1-c)={mp.nstr(Gc/(1-c),12)}")
    for lV in [4, 6, 8]:
        V = mp.mpf(10)**lV
        uV = uV_solve(V)
        lhs = (I(V, c) - (1-c)*uV)*uV
        print(f"   V=1e{lV}: (I-(1-c)u_V)*u_V = {mp.nstr(lhs,10)}   vs G/2? {mp.nstr(Gc,10)}  diff={mp.nstr(lhs-Gc,5)}")
