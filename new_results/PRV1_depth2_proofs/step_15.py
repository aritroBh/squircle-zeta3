
import mpmath as mp
mp.mp.dps = 50

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

def H_dec(t): return depth2(1-t, -t/(1-t))
def E(t):
    return (-mp.polylog(3, -t) + mp.polylog(3, 1-t**2) - 2*mp.polylog(3, 1+t)
            - 2*mp.log(1+t)*mp.polylog(2, 1-t))

for t0 in ['0.2','0.4','0.618','0.8']:
    t = mp.mpf(t0)
    d = H_dec(t) - E(t)
    print('t=', t0, ' H-E =', mp.nstr(d, 30))
