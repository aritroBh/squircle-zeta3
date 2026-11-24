import mpmath as mp
mp.mp.dps = 60
def Li3(x): return mp.re(mp.polylog(3, mp.mpf(x)))
def Li2(x): return mp.re(mp.polylog(2, mp.mpf(x)))
def G3(w):
    """G(a,b,c;1) convergent, via nested quadrature."""
    a,b,c = [mp.mpf(v) for v in w]
    f = lambda t3: mp.quad(lambda t2: mp.quad(lambda t1: 1/(t1-a), [0, t2])/(t2-b), [0, t3])/(t3-c)
    return mp.quad(f, [0,1])
def Nint(x, y):
    """integral int_0^1 log^2((t-x)/(-x))/(t-y) dt, direct."""
    x,y = mp.mpf(x), mp.mpf(y)
    return mp.quad(lambda t: mp.log((t-x)/(-x))**2/(t-y), [0, 1])
def Nform(x, y):
    """closed form: (x/(x-y))*[2 Li3(c) + log^2 X log(1-cX) + 2 log X Li2(cX) - 2 Li3(cX)]"""
    x,y = mp.mpf(x), mp.mpf(y)
    c = x/(x-y); X = (x-1)/x; cX = (x-1)/(x-y)
    lX = mp.log(X) if X>0 else mp.log(abs(X))  # handle sign below via complex
    # use complex logs to be safe
    Xc = mp.mpc(X); 
    lX = mp.log(Xc)
    return (x/(x-y))*(2*mp.polylog(3,c) + lX**2*mp.log(1-cX) + 2*lX*mp.polylog(2,cX) - 2*mp.polylog(3,cX))
for (x,y) in [(2,3),(3,2),(2,-1),(1/2,3),(4,1/3),(2,1/2)]:
    g = G3((x,x,y)); nh = Nint(x,y); nf = Nform(x,y)
    print(f'(x,y)=({x},{y}): G(x,x,y)={mp.nstr(g,12)}  1/2 N direct={mp.nstr(nh/2,12)}  1/2 N closed={mp.nstr(nf/2,12)}  diffs: {mp.nstr(g-nh/2,3)}, {mp.nstr(mp.re(g-nf/2),3)}')
