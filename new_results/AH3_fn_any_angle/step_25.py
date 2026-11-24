
import mpmath as mp
mp.mp.dps = 50

def hc(s, c):
    return 1 + s - mp.sqrt(1 + 2*c*s + s**2)

def I(V, c):
    return mp.quad(lambda u: hc(V*u*mp.exp(-u**2), c), [0, 1, 3, mp.inf], maxdegree=8)

def uV_solve(V):
    return mp.findroot(lambda u: u**2 + mp.log(u) - mp.log(V), 2.0)

for c in [mp.mpf(0), mp.mpf('0.5'), mp.mpf('-0.5')]:
    print(f"c={c}")
    for V in [mp.mpf(20), mp.mpf(100), mp.mpf(702.48), mp.mpf(3000)]:
        uV = uV_solve(V)
        IV = I(V, c)
        # candidates: I(V) = (1-c) u_V + K(c) + o(1)?
        print(f"  V={mp.nstr(V,8)}: I(V)={mp.nstr(IV,12)}  u_V={mp.nstr(uV,10)}  I-(1-c)u_V={mp.nstr(IV-(1-c)*uV,10)}")
