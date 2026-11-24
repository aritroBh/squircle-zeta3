
import mpmath as mp
mp.mp.dps = 60
pi = mp.pi

def psi2(s,h):
    return (1+mp.e**(-s*h/(1-h)))**(-1-h)

def coeff_psi2(s, m, rho=None):
    # Cauchy integral around |h|=rho
    if rho is None:
        rho = 0.8*pi/mp.sqrt(s**2+pi**2)
    f = lambda t: psi2(s, rho*mp.e**(1j*t))/(rho*mp.e**(1j*t))**m
    return mp.quad(f, [0, 2*pi])/(2*pi)

for s in [10.0, 30.0]:
    for m in [9, 10]:  # m+1 even / odd
        cm = coeff_psi2(s,m)
        # asymptotic: (1/s)*2*Re((1+i s/pi)^{m+1})   [pair sum with residue 1/s]
        asym = (1/s)*2*mp.re((1+1j*s/pi)**(m+1))
        print(f"s={s} m={m}: [h^m]psi2 = {mp.nstr(cm,12)}   asym = {mp.nstr(asym,12)}   ratio = {mp.nstr(cm/asym,6)}")
