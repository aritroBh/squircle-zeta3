# step_49: investigate the golden odd part: chi3(1/phi) vs L(chi5,3)
import mpmath as mp
mp.mp.dps = 150
def Li3(x): return mp.re(mp.polylog(3, mp.mpf(x)))
def Li2(x): return mp.re(mp.polylog(2, mp.mpf(x)))
def chi5(n):
    r = n % 5
    return 0 if r==0 else (1 if r in (1,4) else -1)
L5 = mp.nsum(lambda n: chi5(int(n))/mp.mpf(n)**3, [1, mp.inf])
z3 = mp.zeta(3)
phi = (1+mp.sqrt(5))/2
lphi = mp.log(phi)
c3 = (Li3(1/phi)-Li3(-1/phi))/2
print("chi3(1/phi) =", c3)
print("chi3/L5 =", c3/L5)
# pslq chi3(1/phi) against L5, pi^2 lphi, lphi^3, zeta3
rel = mp.pslq([c3, L5, mp.pi**2*lphi, lphi**3, z3], tol=mp.mpf('1e-120'), maxcoeff=10**25)
print("pslq [chi3, L5, pi^2 lphi, lphi^3, z3]:", rel)
# and D = X1-Xm1 against L5 ONLY (2-term, very high precision)
def L3(x):
    x = mp.mpf(x)
    return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(1-x)
D = L3(1/phi) - L3(-1/phi)
print("D =", D)
print("D/L5 =", D/L5)
# D = 2 chi3 + lphi pi^2/6 - lphi^3/2 ?
D2 = 2*c3 + lphi*mp.pi**2/6 - lphi**3/2
print("D - (2chi3 + lphi pi^2/6 - lphi^3/2):", D-D2)
# test classical guess: chi3(1/phi) = 3/10 L5? etc
for q in ['0.7','0.72','0.728','0.73','0.75','0.8']:
    print(q, abs(mp.mpf(q)-c3/L5))
