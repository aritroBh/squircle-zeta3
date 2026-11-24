
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

D = lambda z: mp.polylog(2, z)

# 1) Li_{1,1} closed form:  Li11(x,y) = log(1-x)log(1-y) + Li2(-y/(1-y)) - Li2(-y(1-x)/(1-y))
def Li11_closed(x, y):
    return mp.log(1-x)*mp.log(1-y) + D(-y/(1-y)) - D(-y*(1-x)/(1-y))

for x, y in [(mp.mpf('0.3'), mp.mpf('0.4')), (mp.mpf('-0.3'), mp.mpf('-0.5')), (mp.mpf('0.7'), mp.mpf('-1.5'))]:
    if max(abs(x), abs(x*y)) < 1:
        r = depth2(x, y, 1) - Li11_closed(x, y)
        print('Li11 closed residual', x, y, mp.nstr(r, 10))

# 2) Li_{1,1} along H_dec curve: x=1-t, y=-t/(1-t): claim Li11 = -log t log(1-t) + Li2(t) - Li2(t^2)
for t0 in ['0.2','0.4','0.61','0.8']:
    t = mp.mpf(t0)
    x, y = 1-t, -t/(1-t)
    lhs = depth2(x, y, 1)
    rhs = -mp.log(t)*mp.log(1-t) + D(t) - D(t**2)
    print('Li11 on H curve, t=', t0, ' residual=', mp.nstr(lhs-rhs, 10))

# 3) H_dec'(t): claim  H' = [log t log(1-t) - Li2(t) + Li2(t^2)]/(1-t) - Li2(1-t)/(1-t) - Li2(-t)/t
def H_dec(t): return depth2(1-t, -t/(1-t))
def H_dec_prime(t):
    return (mp.log(t)*mp.log(1-t) - D(t) + D(t**2))/(1-t) - D(1-t)/(1-t) - D(-t)/t
for t0 in ['0.2','0.4','0.61','0.8']:
    t = mp.mpf(t0)
    r = mp.diff(H_dec, t) - H_dec_prime(t)
    print('H_dec derivative, t=', t0, ' residual=', mp.nstr(r, 10))
