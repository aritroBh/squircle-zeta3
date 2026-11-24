import mpmath as mp
mp.mp.dps = 50
m = 2; A = mp.pi/m
def S(th, d):
    p = 1+d
    return mp.fsum([mp.sin(th + i*A)**p for i in range(m)])
def rho_num(th, d):
    p = 1+d
    return S(th, d)**(-1/p)
th = mp.mpf('0.3')
h = mp.mpf('1e-9')
f = lambda dd: rho_num(th, dd)
d1 = (f(-2*h) - 8*f(-h) + 8*f(h) - f(2*h))/(12*h)
sv = [mp.sin(th+i*A) for i in range(m)]
ls = [mp.log(v) for v in sv]
T = [mp.fsum([sv[i]*ls[i]**k for i in range(m)]) for k in range(3)]
S0 = T[0]; U1 = T[1]/S0; U2 = T[2]/(2*S0) - T[1]**2/(2*S0**2)
R1 = mp.log(S0) - U1
R2 = -mp.log(S0) + U1 - U2
print("numeric rho1*rho0:", mp.nstr(d1, 25), " formula:", mp.nstr((1/S0)*R1, 25))
d2 = (2*f(-2*h) - 32*f(-h) + 60*f(0) - 32*f(h) + 2*f(2*h))/(24*h*h)  # wait use standard 2nd deriv stencil
d2 = (-f(-2*h) + 16*f(-h) - 30*f(0) + 16*f(h) - f(2*h))/(12*h*h)
print("numeric rho'' :", mp.nstr(d2, 25), " formula 2*rho0*rho2:", mp.nstr(2*(1/S0)*(R2 + R1**2/2), 25))
