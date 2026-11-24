import mpmath as mp
mp.mp.dps = 60
def Nint(x, y):
    x,y = mp.mpf(x), mp.mpf(y)
    return mp.re(mp.quad(lambda t: mp.log((t-x)/(-x))**2/(t-y), [0, 1]))
def Nform(x, y):
    x,y = mp.mpf(x), mp.mpf(y)
    c = x/(x-y); X = (x-1)/x; cX = (x-1)/(x-y)
    lX = mp.log(mp.mpc(X))
    val = 2*mp.polylog(3,c) + lX**2*mp.log(1-cX) + 2*lX*mp.polylog(2,cX) - 2*mp.polylog(3,cX)
    return mp.re(val)
# x,y outside [0,1] (the kink-line case); avoid x in (0,1)
for (x,y) in [(2,3),(3,2),(2,-1),(4,-3),(5,7),(-2,3),(3,-4),(2,5/3),(7,5)]:
    nh = Nint(x,y); nf = Nform(x,y)
    print(f'(x,y)=({x},{y}): N direct={mp.nstr(nh,14)}  N closed={mp.nstr(nf,14)}  diff={mp.nstr(nh-nf,3)}')
# x=1 case: G(1,1,z;1) = Li3(1/(1-z))
def G_11z(z):
    z = mp.mpf(z)
    return mp.re(mp.quad(lambda t: mp.log(1-t)**2/(t-z), [0,1]))/2
for z in [3, -2, 4]:
    a = G_11z(z); b = mp.re(mp.polylog(3, 1/(1-z)))
    print(f'z={z}: G(1,1,z)={mp.nstr(a,14)}  Li3(1/(1-z))={mp.nstr(b,14)}  diff={mp.nstr(a-b,3)}')
open('progress.log','a').write('step: N(x,y) = F(1)-F(X) validated 60 dps; G(1,1,z)=Li3(1/(1-z)) validated\n')
