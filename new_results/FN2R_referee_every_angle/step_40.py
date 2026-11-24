import mpmath as mp, time
mp.mp.dps = 60
t0=time.time()
def hc(s, c):
    if s <= 0: return mp.mpf(0)
    r = mp.sqrt(1 + 2*c*s + s**2)
    return (1-c) - (1-c**2)/(r + s + c)
f  = lambda x: (1+x)/2*mp.exp(x-1)
fp = lambda x: (2+x)/2*mp.exp(x-1)
def Dn_bump(th, n):
    nn = mp.mpf(n); c = mp.cos(th); lnn = mp.log(nn)
    def bump(v):
        tt = mp.exp(v)/nn
        if tt > 1: return mp.mpf(0)
        fx = f(1-tt)
        if fx <= 0: return mp.mpf(0)
        lyp = mp.log(nn) + (nn-1)*mp.log(fx) + mp.log(fp(1-tt))
        return hc(mp.exp(lyp), c)*mp.exp(v)/nn
    pts = sorted(set([-90, -30, lnn-8, lnn-2, lnn+2, lnn]))
    return mp.quad(bump, pts, maxdegree=10)
def Dn_via_length(th, n):
    nn = mp.mpf(n); c = mp.cos(th); lnn = mp.log(nn)
    def integrand(v):
        tt = mp.exp(v)/nn
        fx = f(1-tt)
        if fx <= 0: return 1*mp.exp(v)/nn
        lyp = mp.log(nn) + (nn-1)*mp.log(fx) + mp.log(fp(1-tt))
        yp = mp.exp(lyp) if lyp < 13000 else mp.exp(13000)
        return mp.sqrt(1+2*c*yp+yp**2)*mp.exp(v)/nn
    pts = sorted(set([-90, -30, lnn-8, lnn-2, lnn]))
    L = mp.quad(integrand, pts, maxdegree=10)
    return 2 - L - f(mp.mpf(0))**n
for thn, th in [('pi/3', mp.pi/3), ('2pi/3', 2*mp.pi/3), ('pi/2', mp.pi/2)]:
    for n in [10**4, 10**6]:
        D1 = Dn_bump(th,n); D2 = Dn_via_length(th,n)
        print("linexp th=%5s n=%7d: |bump - (2-L-f0^n)| = %s" % (thn,n,mp.nstr(abs(D1-D2),3)))
print("elapsed %.1fs"%(time.time()-t0))
