import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); a = mp.pi/(2*sqrt2); pi = mp.pi
bb = mp.e**(1j*mp.pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))

import random
random.seed(41)
w37 = mp.mpf('0'); wh = mp.mpf('0'); w37b = mp.mpf('0')
for _ in range(300):
    y = mp.mpf(10)**random.uniform(-8,-0.4)
    N = pi/(2*y*y); th = pi/(sqrt2*y)
    f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
    rho = mp.findroot(f, 1+y*y/2, tol=mp.mpf('1e-50'))
    sig = (rho-1)/y**2
    x0 = 1/(rho*y); t = rho*y
    phi_true = mp.im(Phi(x0,N)) - th
    zz = sqrt2*rho*y/(2+sqrt2*rho*y); z2 = zz*zz
    A1 = 1 - mp.atan(zz)/zz          # in [0, z2/(1+z2)] as derived
    # step_37 Bh:
    T2n = 2*y*y*sig - sqrt2*rho*y - 2*rho*A1*y
    Bh37 = -y*sig/(2*rho) + T2n/(2*y*(2+sqrt2*rho*y))
    # D13h breg (phi = (pi/rt2)*y*breg):
    breg = -sig/(2*rho) + (-2*y*sig + 2*sqrt2*sig + sqrt2*rho - 2*sqrt2*rho*A1)/(2*sqrt2*(2+sqrt2*rho*y))
    r37 = abs(phi_true - (pi/sqrt2)*Bh37)
    rh_ = abs(phi_true - (pi/sqrt2)*y*breg)
    w37 = max(w37, r37); wh = max(wh, rh_)
    # also verify A1 ∈ [0, z2/(1+z2)]
    w37b = max(w37b, A1 - z2/(1+z2))
print("phi identity: max |phi_true - (pi/rt2) Bh(step_37)| =", mp.nstr(w37,8))
print("phi identity: max |phi_true - (pi/rt2) y breg(D13h)| =", mp.nstr(wh,8))
print("A1 enclosure: max(A1 - z2/(1+z2)) (must be <0):", mp.nstr(w37b,8), " ; A1 >= 0 by arctan<=z")
# analytic: phi_true/(pi/rt2) should equal Bh37; print a sample
y = mp.mpf('0.1'); N = pi/(2*y*y); th = pi/(sqrt2*y)
f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
rho = mp.findroot(f, 1+y*y/2, tol=mp.mpf('1e-50')); sig=(rho-1)/y**2; x0=1/(rho*y); t=rho*y
phi_true = mp.im(Phi(x0,N)) - th
zz = sqrt2*t/(2+sqrt2*t); A1 = 1 - mp.atan(zz)/zz
T2n = 2*y*y*sig - sqrt2*rho*y - 2*rho*A1*y
Bh37 = -y*sig/(2*rho) + T2n/(2*y*(2+sqrt2*t))
print("sample: phi_true =", mp.nstr(phi_true,14), " (pi/rt2)Bh37 =", mp.nstr((pi/sqrt2)*Bh37,14))
