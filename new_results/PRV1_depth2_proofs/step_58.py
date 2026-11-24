import mpmath as mp
mp.mp.dps = 150

def depth2(x, y, p=2):
    rho = max(abs(x), abs(x*y))
    target = mp.power(10, -mp.mp.dps - 5)
    nmax = int(mp.ceil(mp.log(target*(1-rho))/mp.log(rho)))
    h, xy_pow, tot = mp.mpf(0), x*y, mp.mpf(0)
    for n in range(1, nmax+1):
        tot += h/n**p
        h = x*h + x*xy_pow/n
        xy_pow *= x*y
    return tot

V = depth2(mp.mpf(-1)/2, mp.mpf(-2)/3)   # c=2/3 family value at t=1/2
print("Li21(-1/2,-2/3) =", mp.nstr(V, 30))

li3 = lambda z: mp.polylog(3, z)
l2, l3 = mp.log(2), mp.log(3)
z2, z3 = mp.pi**2/6, mp.zeta(3)
basis = [z3, li3(mp.mpf(1)/2), li3(mp.mpf(-1)/2), li3(mp.mpf(1)/3), li3(mp.mpf(-1)/3),
         li3(mp.mpf(2)/3), li3(mp.mpf(-2)/3), li3(mp.mpf(2)/5), li3(mp.mpf(-2)/5),
         li3(mp.mpf(3)/5), li3(mp.mpf(-3)/5), li3(mp.mpf(1)/4), li3(mp.mpf(3)/4),
         z2*l2, z2*l3, l2**3, l3**3, l2**2*l3, l2*l3**2]
names = ['z3','L3(1/2)','L3(-1/2)','L3(1/3)','L3(-1/3)','L3(2/3)','L3(-2/3)',
         'L3(2/5)','L3(-2/5)','L3(3/5)','L3(-3/5)','L3(1/4)','L3(3/4)',
         'z2l2','z2l3','l2^3','l3^3','l2^2l3','l2l3^2']
rel = mp.pslq([V]+basis, tol=mp.power(10,-120), maxcoeff=mp.mpf(10)**8, maxsteps=10000)
if rel is None:
    print("PSIQ: NO relation found (basis of 19, 150 dps)")
else:
    print("relation:", rel)
    s = rel[0]*V + sum(r*b for r,b in zip(rel[1:], basis))
    print("residual:", mp.nstr(s, 5))
