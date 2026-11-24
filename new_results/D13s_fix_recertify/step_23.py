
import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); pi = mp.pi; a = pi/(2*sqrt2); b = mp.e**(1j*pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(b*x/(1+b*x))

t = mp.mpf('1e-2')
rho = mp.sqrt(1 + t*t/(1+sqrt2*t))
y = t/rho; N = pi/(2*y*y); th = pi/(sqrt2*y); x0 = 1/(rho*y)
phi = mp.im(Phi(x0,N)) - th
Btrue = (phi + pi/4)/((pi/sqrt2)*y)
print("Btrue =", mp.nstr(Btrue, 14))
print("phi+pi/4 =", mp.nstr(phi+pi/4, 14), " (pi/rt2)y =", mp.nstr((pi/sqrt2)*y, 14))
# term-by-term: phi = a(1-2rho)/t + N*(pi/4 - atan(1/(1+sqrt2 t)))
T1 = a*(1-2*rho)/t
T2 = N*(pi/4 - mp.atan(1/(1+sqrt2*t)))
print("T1 =", mp.nstr(T1,14), " T2 =", mp.nstr(T2,14), " T1+T2 =", mp.nstr(T1+T2,14), " phi =", mp.nstr(phi,14))
# so phi + pi/4 = T1 + T2 + pi/4 ; divide by (pi/sqrt2)y = 2a t/rho:
d = 2*a*t/rho
print("T1/d =", mp.nstr(T1/d,14))
print("T2/d =", mp.nstr(T2/d,14))
print("(pi/4)/d =", mp.nstr((pi/4)/d,14))
print("sum =", mp.nstr(T1/d+T2/d+(pi/4)/d, 14))
