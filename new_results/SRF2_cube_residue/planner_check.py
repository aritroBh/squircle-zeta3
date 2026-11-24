# Planner check of SRF2: rho_{-2} = 3(-2 + I2'), V_40 vs R, and a closed-form attempt for I2'.
import json
from mpmath import mp, mpf, quad, sqrt, gamma, pi, pslq, nstr, log
mp.dps = 40
I2 = quad(lambda a: (1 + a)*a**3/(sqrt(1 + a**6) + 1), [0, 1])   # = int ((1+a)sqrt(1+a^6) - 1 - a)/a^3, cancellation-free
rho2 = 3*(-2 + I2)
R = -1 - 2*pi**mpf(1.5)/gamma(mpf(1)/4)**2
print("I2' =", nstr(I2, 25), " rho_-2 =", nstr(rho2, 20))
V = [mpf(x) for x in json.load(open('EV_coeffs_K40.json'))['V']]
for k in (12, 14, 16):
    print("k =", k, " (V_k - R(-1)^k)(-1)^k 2^(k+1) =", nstr((V[k] - R*(-1)**k)*(-1)**k*2**(k+1), 8))
print("(-1)^40 V_40 - R =", nstr(V[40] - R, 3))
g13, g16 = gamma(mpf(1)/3), gamma(mpf(1)/6)
cands = [I2, 1, sqrt(2), sqrt(3), g13**3/pi, pi**2/g13**3, sqrt(3)*g13**3/pi, sqrt(3)*pi**2/g13**3, log(1+sqrt(2)), sqrt(2)*log(1+sqrt(2))]
print("PSLQ I2' on {1, sqrt2, sqrt3, G(1/3)^3/pi, pi^2/G(1/3)^3, (x sqrt3), log(1+sqrt2), sqrt2 log(1+sqrt2)}:", pslq(cands, maxcoeff=10**5, maxsteps=10**6))
