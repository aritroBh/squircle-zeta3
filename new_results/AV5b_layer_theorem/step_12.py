import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
step = r'''
import mpmath as mp
mp.mp.dps = 60

def make_speed(m):
    def S(phi, p):
        return mp.fsum(abs(mp.cos(phi - i*mp.pi/m))**p for i in range(1, m+1))
    def speed(phi, p):
        Sp = S(phi, p)
        rho = Sp**(-mp.mpf(1)/p)
        # rho' = d rho/d phi
        drho = mp.diff(lambda ph: S(ph, p)**(-mp.mpf(1)/p), phi)
        return mp.sqrt(rho**2 + drho**2)
    return speed

def P_of_p(m, p):
    speed = make_speed(m)
    # integrate over [0, 2*pi), splitting at vertex/edge-mid directions for smoothness
    pts = [mp.pi*k/m for k in range(2*m)] + [mp.pi/2 + mp.pi*k/m for k in range(2*m)]
    pts = sorted(set(mp.nstr(x, 30) for x in pts))
    pts = [mp.mpf(x) for x in pts]
    pts = [x for x in pts if 0 < x < 2*mp.pi]
    total = mp.mpf('0')
    bounds = [mp.mpf('0')] + pts + [2*mp.pi]
    for a, b in zip(bounds[:-1], bounds[1:]):
        total += mp.quad(lambda ph: speed(ph, p), [a, b])
    return total

# p-derivatives of P at p=1 via differentiating the integrand in p
def P_deriv(m, k):
    speed = make_speed(m)
    pts = sorted([mp.pi*kk/m for kk in range(2*m)] + [mp.pi/2 + mp.pi*kk/m for kk in range(2*m)])
    pts = [x for x in pts if 0 < x < 2*mp.pi]
    bounds = [mp.mpf('0')] + pts + [2*mp.pi]
    total = mp.mpf('0')
    for a, b in zip(bounds[:-1], bounds[1:]):
        total += mp.quad(lambda ph: mp.diff(lambda pp: speed(ph, pp), mp.mpf(1), k), [a, b])
    return total

for m in [2, 3]:
    P0 = P_of_p(m, mp.mpf(1))
    print(f"m={m}: P(1) = {mp.nstr(P0, 30)}  (formula 4m(sec-cos) = {mp.nstr(4*m*(1/mp.cos(mp.pi/(2*m))-mp.cos(mp.pi/(2*m))), 30)})")
for m in [2, 3]:
    for k in [1, 2, 3]:
        v = P_deriv(m, k)
        print(f"m={m}: P^({k})(1) = {mp.nstr(v, 25)}")
    print()
# known: m=2 P'''(1) = 6 sqrt2 zeta3 - 2 sqrt2 pi^2
print("m=2 target P'''(1):", mp.nstr(6*mp.sqrt(2)*mp.zeta(3) - 2*mp.sqrt(2)*mp.pi**2, 25))
print("m=3 C_3 zeta(3):", mp.nstr(-9*mp.sqrt(3)/8*mp.zeta(3), 25))
'''
with open("step_05.py","w") as f: f.write(step)
exec(step)
