import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); a = mp.pi/(2*sqrt2); pi = mp.pi
bb = mp.e**(1j*mp.pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))
import random
random.seed(43)
w = mp.mpf('0'); w2 = mp.mpf('0'); wA = mp.mpf('-9')
for _ in range(300):
    y = mp.mpf(10)**random.uniform(-8,-0.4)
    N = pi/(2*y*y); th = pi/(sqrt2*y)
    f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
    rho = mp.findroot(f, 1+y*y/2, tol=mp.mpf('1e-50'))
    sig = (rho-1)/y**2; x0 = 1/(rho*y); t = rho*y
    phi_true = mp.im(Phi(x0,N)) - th
    zz = sqrt2*t/(2+sqrt2*t); z2 = zz*zz
    Ap = 1 - mp.atan(zz)/zz
    A1 = Ap/y                      # the value that makes Bh exact
    T2n = 2*y*y*sig - sqrt2*rho*y - 2*rho*A1*y
    Bh = -y*sig/(2*rho) + T2n/(2*y*(2+sqrt2*t))
    w = max(w, abs(phi_true - (pi/sqrt2)*Bh))
    # D13h breg with A1' (their convention, since A1mreg encloses A1'/y? check both)
    breg_h = -sig/(2*rho) + (-2*y*sig + 2*sqrt2*sig + sqrt2*rho - 2*sqrt2*rho*Ap)/(2*sqrt2*(2+sqrt2*t))
    w2 = max(w2, abs(phi_true - (pi/sqrt2)*y*breg_h))
    wA = max(wA, Ap - z2/(1+z2))   # A1' <= z2/(1+z2) ?
print("phi = (pi/rt2) Bh[A1'/y]  (step_37 formula): worst residual =", mp.nstr(w,10))
print("phi = (pi/rt2) y breg_h[A1'] (D13h formula): worst residual =", mp.nstr(w2,10))
print("A1' <= z2/(1+z2): max(A1' - z2/(1+z2)) =", mp.nstr(wA,8))
# so: step_37's Bh formula is EXACT with A1 = A1'/y ; D13h's breg is EXACT with A1 = A1' ?
# => D13h A1mreg = [0, z2/((1+z2)y)] must contain A1'/y <=> A1' <= z2/(1+z2) (holds)
# => step_37/38 A1m = [0, z2/(1+z2)] must contain A1'/y <=> A1' <= y z2/(1+z2) -- FAILS for y<1/3
print()
print("check: A1'/y vs step_38 A1m upper z2/(1+z2): A1'/y <= z2/(1+z2) ?")
m = mp.mpf('-9')
for _ in range(500):
    y = mp.mpf(10)**random.uniform(-8,-0.4)
    t = y*(1+y*y*random.random())
    zz = sqrt2*t/(2+sqrt2*t); z2=zz*zz
    Ap = 1-mp.atan(zz)/zz
    m = max(m, Ap/y - z2/(1+z2))
print("max(A1'/y - z2/(1+z2)) over random (must be <=0 for step_38 box to be correct):", mp.nstr(m,8), "=> POSITIVE => step_38 A1m box CONFIRMED WRONG")
