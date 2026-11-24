
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

def H_dec(t): return depth2(1-t, -t/(1-t))
def H_dec_closed(t):
    c = (1-t)/(1+t)
    return (-mp.polylog(3, -t) + 2*mp.polylog(3, 1-t) + mp.polylog(3, 1-t**2)/2
            - mp.polylog(3, c) + mp.polylog(3, -c) - 2*mp.log(1+t)*mp.polylog(2, 1-t)
            - 3*mp.zeta(3)/4)

print("H_dec closed form check:")
for t0 in ['0.1','0.2','0.33','0.5','0.61803398874989484820458683436563811772030917980576','0.8','0.9']:
    t = mp.mpf(t0)
    r = H_dec(t) - H_dec_closed(t)
    print('  t=', t0, ' residual=', mp.nstr(r, 8))

# Specialisation at t = 1/phi
phi = (1+mp.sqrt(5))/2
t = 1/phi
l = mp.log(phi); z2 = mp.zeta(2); z3 = mp.zeta(3)
L3p = mp.polylog(3, 1/phi); L3m = mp.polylog(3, -1/phi)
mine = H_dec_closed(t)
dec1 = 3*z3/5 - 7*l*z2/5 + 2*l**3 + L3m + L3p
print("swap_dec: my closed form =", mp.nstr(mine, 40))
print("swap_dec: DEC1 LLL form =", mp.nstr(dec1, 40))
print("difference =", mp.nstr(mine - dec1, 12))
# Also direct series:
print("direct series =", mp.nstr(H_dec(t), 40))
# G_dec at t=1/phi:
def G_dec_closed(t):
    return (z3 - mp.log(t)*mp.polylog(2, t**2)/2
            + mp.polylog(3, t**2) - 2*mp.polylog(3, t) - mp.polylog(3, 1-t))
g = G_dec_closed(t)
lam_dec1 = 2*z3/5 + 4*l*z2/5 - l**3 + 3*L3m + L3p
print("Lambda_dec: my closed form =", mp.nstr(g, 40))
print("Lambda_dec: DEC1 form     =", mp.nstr(lam_dec1, 40))
print("difference =", mp.nstr(g - lam_dec1, 12))
