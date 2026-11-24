
import mpmath as mp
mp.mp.dps = 50

def hc(s, c):
    return 1 + s - mp.sqrt(1 + 2*c*s + s**2)

def I(V, c):
    return mp.quad(lambda u: hc(V*u*mp.exp(-u**2), c), [0, 1, 3, mp.inf], maxdegree=9)

def uV_solve(V):
    return mp.findroot(lambda u: u**2 + mp.log(u) - mp.log(V), (mp.log(V))**mp.mpf('0.5'))

for c in [mp.mpf(0), mp.mpf('0.5'), mp.mpf('-0.5')]:
    print(f"c={c}")
    for lV in [2,3,4,5,6,8]:
        V = mp.mpf(10)**lV
        uV = uV_solve(V)
        IV = I(V, c)
        R = IV/((1-c)*uV) - 1
        print(f"  V=1e{lV}: u_V={mp.nstr(uV,8)} I/(1-c)u_V-1={mp.nstr(R,6)}  R*u_V^2={mp.nstr(R*uV**2,6)} R*u_V={mp.nstr(R*uV,6)}")
