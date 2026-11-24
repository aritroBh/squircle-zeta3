import mpmath as mp, time
mp.mp.dps = 50
t0 = time.time()

def hc(s, c):
    if s <= 0: return mp.mpf(0)
    if s < 1: return 1 + s - mp.sqrt(1 + 2*c*s + s**2)
    r = mp.sqrt(1 + 2*c*s + s**2)
    return (1-c) - (1-c**2)/(r + s + c)
def S12(z):
    return mp.zeta(3) - mp.polylog(3, 1-z) + mp.log(1-z)*mp.polylog(2, 1-z) \
           + mp.mpf('0.5')*mp.log(1-z)**2*mp.log(z)
def cumulants(th):
    z = mp.sin(th/2)**2; w = 1/mp.tan(th/2)**2; ell = mp.log(mp.cos(th/2)**2); u = w*ell
    return (1+u, 1+mp.pi**2/3-2*mp.polylog(2,z)-u**2, -12*w*S12(z)+2*(1+u**3))

def I_model(V, c):
    return mp.quad(lambda u: hc(V*u*mp.exp(-u**2), c), [0, 1, 3, mp.inf], maxdegree=10)
def ustar(V):
    return mp.findroot(lambda u: u**2 - mp.log(u) - mp.log(V), mp.sqrt(mp.log(V)))

print("(I-(1-c)u*)u*  vs  (1-c)k1/2  and  vs  (1-c)k1/2 - (1-c^2)/4")
for cval, thname in [(mp.mpf(0),'pi/2'), (mp.mpf('0.5'),'pi/3'), (mp.mpf('-0.5'),'2pi/3')]:
    th = mp.acos(cval)
    k1 = cumulants(th)[0]
    candA = (1-cval)*k1/2
    candB = candA - (1-cval**2)/4
    print(f"theta={thname}: candA=(1-c)k1/2={mp.nstr(candA,10)}, candB={mp.nstr(candB,10)}")
    for lV in [3, 5, 8, 10]:
        V = mp.mpf(10)**lV
        us = ustar(V)
        lhs = (I_model(V, cval) - (1-cval)*us)*us
        print(f"   V=1e{lV}: (I-(1-c)u*)u*={mp.nstr(lhs,10)}  -candA={mp.nstr(lhs-candA,4)}  -candB={mp.nstr(lhs-candB,4)}")
print("elapsed", time.time()-t0)
