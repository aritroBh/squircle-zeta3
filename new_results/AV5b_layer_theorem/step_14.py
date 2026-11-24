
import mpmath as mp
mp.mp.dps = 60

def make_speed(m):
    def S(phi, p):
        return mp.fsum(abs(mp.cos(phi - i*mp.pi/m))**p for i in range(1, m+1))
    def speed(phi, p):
        Sp = S(phi, p)
        rho = Sp**(-mp.mpf(1)/p)
        drho = mp.diff(lambda ph: S(ph, p)**(-mp.mpf(1)/p), phi)
        return mp.sqrt(rho**2 + drho**2)
    return speed

def sector_bounds(m):
    pts = sorted(set([mp.pi*kk/m for kk in range(2*m)] + [mp.pi/2 + mp.pi*kk/m for kk in range(2*m)]))
    pts = [x for x in pts if 0 < x < 2*mp.pi]
    return [mp.mpf(0)] + pts + [2*mp.pi]

def P_deriv(m, k):
    speed = make_speed(m)
    total = mp.mpf(0)
    for a, b in zip(sector_bounds(m)[:-1], sector_bounds(m)[1:]):
        total += mp.quad(lambda ph: mp.diff(lambda pp: speed(ph, pp), mp.mpf(1), k), [a, b])
    return total

tg2 = 6*mp.sqrt(2)*mp.zeta(3) - 2*mp.sqrt(2)*mp.pi**2
print("m=2 target P3(1):", mp.nstr(tg2, 25))
for m in [2, 3]:
    vals = {}
    for k in [1, 2, 3]:
        vals[k] = P_deriv(m, k)
        print(f"m={m}: P^({k})(1) = {mp.nstr(vals[k], 25)}")
    if m == 2:
        print("  residual vs target:", mp.nstr(vals[3]-tg2, 8))
    if m == 3:
        C3 = -9*mp.sqrt(3)/8
        print("  P3(1)/(zeta(3)) =", mp.nstr(vals[3]/mp.zeta(3), 20))
        print("  C_3 =", mp.nstr(C3, 20), " residual:", mp.nstr(vals[3]-C3*mp.zeta(3), 8))
        print("  P1(1) target 2 sqrt3 ln2:", mp.nstr(2*mp.sqrt(3)*mp.log(2), 25))
