
import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); pi = mp.pi; a = pi/(2*sqrt2); b = mp.e**(1j*pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(b*x/(1+b*x))

# Direct Btrue at true saddles for small t
for tv in ['1e-1','1e-2','1e-3','1e-4','1e-6']:
    t = mp.mpf(tv)
    rho = mp.sqrt(1 + t*t/(1+sqrt2*t))
    y = t/rho; N = pi/(2*y*y); th = pi/(sqrt2*y); x0 = 1/(rho*y)
    phi = mp.im(Phi(x0,N)) - th
    Btrue = (phi + pi/4)/((pi/sqrt2)*y)
    # my derived closed form:
    Bform = rho*(1-rho)**2/(2*t**2) + rho**3*(pi/4 - mp.atan(1/(1+sqrt2*t)))/(sqrt2*t**3) - rho**3/(2*t**2) + sqrt2*rho/(4*t) - rho*(1-rho)**2/(2*t**2)
    print(f"t={tv}: Btrue={mp.nstr(Btrue,14)}")

# derive carefully: phi = a(1-2rho)/t + N[pi/4 - atan(1/(1+sqrt2 t))]
print()
print("check phi decomposition:")
for tv in ['1e-2','1e-4']:
    t = mp.mpf(tv)
    rho = mp.sqrt(1 + t*t/(1+sqrt2*t))
    y = t/rho; N = pi/(2*y*y); th = pi/(sqrt2*y); x0 = 1/(rho*y)
    phi = mp.im(Phi(x0,N)) - th
    phi2 = a*(1-2*rho)/t + N*(pi/4 - mp.atan(1/(1+sqrt2*t)))
    print(f"t={tv}: phi={mp.nstr(phi,14)}  phi2={mp.nstr(phi2,14)}  resid={mp.nstr(abs(phi-phi2),6)}")
