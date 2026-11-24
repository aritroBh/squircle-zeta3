import mpmath as mp
mp.mp.dps = 60
def phi(s,c):
    if s==0: return 1-c
    r = mp.sqrt(1+2*c*s+s**2)
    return (1-c**2)/(r*(r+s+c))
th = mp.pi/6; c = mp.cos(th); zv = mp.sin(th/2)**2
for wv in [mp.mpf('-0.9'), mp.mpf('0.9')]:
    # s = u^p smoothing for the s in [0,1] part; p chosen so (u^p)^w u^{p-1} ~ u^{p(w+1)-1} smooth
    p = 10
    i1 = mp.quad(lambda u: (u**p)**wv*phi(u**p,c)*p*u**(p-1), [0, 1])
    # tail: s in [1,inf), s = u^{-p}: s^w = u^{-pw}, ds = -p u^{-p-1}
    i2 = mp.quad(lambda u: (u**-p)**wv*phi(u**-p,c)*p*u**(-p-1), [0, 1])
    lhs = (i1+i2)/(1-c)
    rhs = (mp.pi*wv/mp.sin(mp.pi*wv))*2*zv*mp.hyper([-wv,1+wv],[2],zv)/(1-c)
    print("w=", wv, " smoothed diff:", mp.nstr(abs(lhs-rhs),3))
