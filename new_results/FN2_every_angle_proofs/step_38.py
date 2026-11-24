import mpmath as mp, time
mp.mp.dps = 70
t0 = time.time()
def hc(s, c):
    if s <= 0: return mp.mpf(0)
    if s < 1: return 1 + s - mp.sqrt(1 + 2*c*s + s**2)
    r = mp.sqrt(1 + 2*c*s + s**2)
    return (1-c) - (1-c**2)/(r + s + c)
def make_f(which):
    if which=='exp':
        e = mp.e; a = e/(e-1)
        return (lambda x: (mp.exp(x)-1)/(e-1), lambda x: mp.exp(x)/(e-1), (a,a,a))
    if which=='quad':
        return (lambda x: (x+x**2)/2, lambda x: (1+2*x)/2, (mp.mpf('1.5'), mp.mpf(1), mp.mpf(0)))

def Dn_cfg(which, th, n, pts, maxdeg=10):
    f, fp, _ = make_f(which)
    nn = mp.mpf(n); c = mp.cos(th)
    def integrand(v):
        t = mp.exp(v)/nn
        if t >= 1: return mp.mpf(0)
        fx = f(1-t)
        if fx <= 0: return mp.mpf(0)
        yp = nn*mp.exp((nn-1)*mp.log(fx))*fp(1-t)
        return hc(yp, c)*mp.exp(v)/nn
    return mp.quad(integrand, pts, maxdegree=maxdeg, error=True)

th = 2*mp.pi/3; n = 10**6
mid = mp.log(mp.log(mp.mpf(n)))
configs = {
 'A [-160,mid-15,mid+10,lnn]': [-160, mid-15, min(mid+10, mp.log(n)-mp.mpf('0.5')), mp.log(n)],
 'B [-160,mid-8,mid+8,lnn]':   [-160, mid-8, min(mid+8, mp.log(n)-mp.mpf('0.5')), mp.log(n)],
 'C [-60,mid-12,mid+12,lnn]':  [-60, mid-12, min(mid+12, mp.log(n)-mp.mpf('0.5')), mp.log(n)],
 'D narrow peak':               [-160, mid-6, mid-2, mid+2, mid+6, mp.log(n)],
}
for name, pts in configs.items():
    v, err = Dn_cfg('exp', th, n, pts)
    print(f"{name}: {mp.nstr(v,30)}  err~{mp.nstr(err,2)}")
print("elapsed", time.time()-t0)
