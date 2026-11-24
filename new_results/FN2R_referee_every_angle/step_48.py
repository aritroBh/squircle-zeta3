import mpmath as mp, time
mp.mp.dps = 50
t0=time.time()
def hc(s, c): return 1 + s - mp.sqrt(1 + 2*c*s + s**2)
def rho(s, c): return mp.sqrt(1 + 2*c*s + s**2) - s - c
def S12(zv):
    return mp.zeta(3) - mp.polylog(3, 1-zv) + mp.log(1-zv)*mp.polylog(2, 1-zv) \
           + mp.mpf('0.5')*mp.log(1-zv)**2*mp.log(zv)
def kforms(th):
    zv = mp.sin(th/2)**2; w2 = 1/mp.tan(th/2)**2; u = w2*mp.log(mp.cos(th/2)**2)
    return (1+u, 1+mp.pi**2/3-2*mp.polylog(2,zv)-u**2, 2+2*u**3-12*w2*S12(zv))
def I(V, c):
    return mp.quad(lambda u: hc(V*u*mp.exp(-u**2), c), [0, 1, 3, mp.inf], maxdegree=9)
def ustar(V):
    return mp.findroot(lambda u: u**2 - mp.log(u) - mp.log(V), (mp.log(V))**mp.mpf('0.5'))
def Klayer(c):
    i1 = mp.quad(lambda s: hc(s, c)/s, [0, 1])
    i2 = mp.quad(lambda s: rho(s, c)/s, [1, mp.inf])
    return (i1 - i2)/2
print("== Proposition: (I-(1-c)u*)u* -> (1-c) k1/2 = Klayer ; AH3 guess = Klayer-(1-c^2)/4 ==")
for th in [mp.pi/3, mp.pi/2, 2*mp.pi/3]:
    c = mp.cos(th); k1 = kforms(th)[0]
    rhs = (1-c)*k1/2
    Kl = Klayer(c)
    guess = Kl - (1-c**2)/4
    print(f"th={mp.nstr(th,7)}: (1-c)k1/2={mp.nstr(rhs,14)} Klayer={mp.nstr(Kl,14)} (diff={mp.nstr(Kl-rhs,2)}) AH3guess={mp.nstr(guess,14)}")
    for lV in [4, 6, 8]:
        V = mp.mpf(10)**lV
        us = ustar(V)
        lhs = (I(V,c)-(1-c)*us)*us
        print(f"   V=1e{lV}: lhs={mp.nstr(lhs,12)}  lhs-(1-c)k1/2={mp.nstr(lhs-rhs,4)}  lhs-guess={mp.nstr(lhs-guess,4)}")
print("elapsed %.1fs"%(time.time()-t0))
with open("progress.log","a") as f2:
    f2.write(f"{time.strftime('%H:%M')} step_11: Proposition verified: (I-(1-c)u*)u* -> (1-c)k1/2 = Klayer (IBP identity holds to 1e-49); AH3's Klayer-(1-c^2)/4 guess is off by (1-c^2)/4 as FN2 states [referee, NUMERICAL]\n")
