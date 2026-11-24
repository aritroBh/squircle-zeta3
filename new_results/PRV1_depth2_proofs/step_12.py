
import mpmath as mp
mp.mp.dps = 60

def depth2(x, y, p=2):
    rho = max(abs(x), abs(x*y))
    assert 0 < rho < 1
    target = mp.power(10, -mp.mp.dps - 5)
    nmax = int(mp.ceil(mp.log(target*(1-rho))/mp.log(rho)))
    h, xy_power, total = mp.mpf(0), x*y, mp.mpf(0)
    for n in range(1, nmax+1):
        total += h/n**p
        h = x*h + x*xy_power/n
        xy_power *= x*y
    return total

def G_dec(t): return depth2(-t, -(1-t)/t)
def G_dec_closed(t):
    return (mp.zeta(3) - mp.log(t)*mp.polylog(2, t**2)/2
            + mp.polylog(3, t**2) - 2*mp.polylog(3, t) - mp.polylog(3, 1-t))

for t0 in ['0.2','0.31','0.5','0.618033988749894848204586834365638117720309179805762862135448','0.8','0.95']:
    t = mp.mpf(t0)
    r = G_dec(t) - G_dec_closed(t)
    print('t=', mp.nstr(t,6), ' residual=', mp.nstr(r, 12))
