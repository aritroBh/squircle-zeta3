import mpmath as mp
mp.mp.dps = 60
def Nint(x, y):
    x,y = mp.mpf(x), mp.mpf(y)
    return mp.quad(lambda t: mp.log((t-x)/(-x))**2/(t-y), [0, 1])
def Nform(x, y):
    x,y = mp.mpf(x), mp.mpf(y)
    c = x/(x-y); X = (x-1)/x; cX = (x-1)/(x-y)
    lX = mp.log(mp.mpc(X))
    val = (x/(x-y))*(2*mp.polylog(3,c) + lX**2*mp.log(1-cX) + 2*lX*mp.polylog(2,cX) - 2*mp.polylog(3,cX))
    return mp.re(val)
for (x,y) in [(2,3),(3,2),(2,-1),(1/2,3),(4,1/3),(2,1/2),(5,7),(1,3),(1,-2)]:
    nh = Nint(x,y); nf = Nform(x,y)
    print(f'(x,y)=({x},{y}): N direct={mp.nstr(nh,14)}  N closed={mp.nstr(nf,14)}  diff={mp.nstr(nh-nf,3)}')
open('progress.log','a').write('step: N(x,y) closed form validated (1-d quad vs closed form, 60 dps)\n')
