import mpmath as mp
mp.mp.dps = 250
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
V = depth2(mp.mpf(-1)/2, mp.mpf(-2)/3)
li3 = lambda z: mp.polylog(3, z)
l2, l3 = mp.log(2), mp.log(3)
z2, z3 = mp.pi**2/6, mp.zeta(3)
# surviving basis after drops: ['z3','L3(1/2)','L3(-1/2)','L3(1/3)','L3(-1/3)','L3(2/3)','L3(-2/3)',
#   'L3(2/5)','L3(-2/5)','L3(3/5)','L3(-3/5)','L3(3/4)','z2l2','z2l3']
basis = [z3, li3(mp.mpf(1)/2), li3(mp.mpf(-1)/2), li3(mp.mpf(1)/3), li3(mp.mpf(-1)/3),
         li3(mp.mpf(2)/3), li3(mp.mpf(-2)/3), li3(mp.mpf(2)/5), li3(mp.mpf(-2)/5),
         li3(mp.mpf(3)/5), li3(mp.mpf(-3)/5), li3(mp.mpf(3)/4), z2*l2, z2*l3]
rel = [51657984, 8819907, -2309587, -1259719, -44065700, -29158608, -4194165,
       -59973542, -37298437, 19825195, -48559972, 5705039, -4545108, 5160575, 11216678]
s = rel[0]*V + sum(r*b for r,b in zip(rel[1:], basis))
print("same combination at 250 dps:", mp.nstr(s, 10))
print("=> spurious" if abs(s) > mp.power(10,-150) else "=> genuine")
