import mpmath as mp
mp.mp.dps = 60
def phi(s,c):
    if s==0: return 1-c
    r = mp.sqrt(1+2*c*s+s**2)
    return (1-c**2)/(r*(r+s+c))
print("per (theta,w) residuals, split quadrature:")
worst=(0,None)
for th in [mp.pi/6, mp.pi/3, mp.pi/2, 2*mp.pi/3, 5*mp.pi/6]:
    c = mp.cos(th); zv = mp.sin(th/2)**2
    for wv in [mp.mpf('-0.9'), mp.mpf('-0.5'), mp.mpf('-0.1'), mp.mpf('0.1'), mp.mpf('0.5'), mp.mpf('0.9')]:
        i1 = mp.quad(lambda s: s**wv*phi(s,c), [0, 1])
        # tail with substitution s=1/t for accuracy
        i2 = mp.quad(lambda t: t**(-wv-2)*phi(1/t,c), [0, 1])
        lhs = (i1+i2)/(1-c)
        rhs = (mp.pi*wv/mp.sin(mp.pi*wv))*2*zv*mp.hyper([-wv,1+wv],[2],zv)/(1-c)
        d = abs(lhs-rhs)
        if d>worst[0]: worst=(d,(th,wv,lhs,rhs))
        print(f" th={mp.nstr(th,7)} w={wv:+}: diff={mp.nstr(d,2)}")
print("worst:", worst[0], worst[1][0], worst[1][1], mp.nstr(worst[1][2],20), mp.nstr(worst[1][3],20))
