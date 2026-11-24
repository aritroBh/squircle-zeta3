
import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); pi = mp.pi; a = pi/(2*sqrt2); b = mp.e**(1j*pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(b*x/(1+b*x))

# single concrete point, arbitrary rho (like review step_31)
y = mp.mpf('0.037'); rho = mp.mpf('1.0011'); sig = (rho-1)/y**2
N = pi/(2*y*y); th = pi/(sqrt2*y); x0 = 1/(rho*y); t = rho*y
phi_true = mp.im(Phi(x0,N)) - th
zz = sqrt2*t/(2+sqrt2*t)
Ap = 1 - mp.atan(zz)/zz
breg_h = -sig/(2*rho) + (-2*y*sig + 2*sqrt2*sig + sqrt2*rho - 2*sqrt2*rho*Ap)/(2*sqrt2*(2+sqrt2*t))
Btrue = (phi_true + pi/4)/((pi/sqrt2)*y)
print("Btrue =", mp.nstr(Btrue, 20))
print("breg_h =", mp.nstr(breg_h, 20))
print("Btrue - breg_h =", mp.nstr(Btrue - breg_h, 20))
print("Btrue - breg_h + 1/12 =", mp.nstr(Btrue - breg_h + mp.mpf(1)/12, 20))
