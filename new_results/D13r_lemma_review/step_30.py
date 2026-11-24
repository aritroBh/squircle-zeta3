import mpmath as mp, random
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); a = mp.pi/(2*sqrt2); pi = mp.pi
bb = mp.e**(1j*mp.pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))
random.seed(44)
w = mp.mpf('0'); wbox = mp.mpf('-9')
for _ in range(400):
    y = mp.mpf(10)**random.uniform(-8,-0.4)
    N = pi/(2*y*y); th = pi/(sqrt2*y)
    f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
    rho = mp.findroot(f, 1+y*y/2, tol=mp.mpf('1e-50'))
    sig = (rho-1)/y**2; x0 = 1/(rho*y); t = rho*y
    phi_true = mp.im(Phi(x0,N)) - th
    zz = sqrt2*t/(2+sqrt2*t); z2=zz*zz
    Ap = 1 - mp.atan(zz)/zz
    breg_h = -sig/(2*rho) + (-2*y*sig + 2*sqrt2*sig + sqrt2*rho - 2*sqrt2*rho*Ap)/(2*sqrt2*(2+sqrt2*t))
    Btrue = (phi_true + pi/4)/((pi/sqrt2)*y)
    w = max(w, abs(Btrue - breg_h))
    # box: A1mreg = [0, 2R^2 Y/((2+rt2 R Y)^2 (1+z2))]; breg affine in A1 with slope -rho/(2+rt2 t) <0
    Aup = 2*rho**2*y/((2+sqrt2*t)**2*(1+z2))
    blo = -sig/(2*rho) + (-2*y*sig + 2*sqrt2*sig + sqrt2*rho - 2*sqrt2*rho*Aup)/(2*sqrt2*(2+sqrt2*t))
    bhi = -sig/(2*rho) + (-2*y*sig + 2*sqrt2*sig + sqrt2*rho)/(2*sqrt2*(2+sqrt2*t))
    ok = blo <= Btrue <= bhi
    if not ok:
        viol = max(blo-Btrue, Btrue-bhi)
        wbox = max(wbox, viol)
print("D13h: (phi+pi/4)/((pi/rt2) y) = breg_h[A1']: worst residual =", mp.nstr(w,10))
print("D13h: Btrue inside breg-box over all samples: max violation =", mp.nstr(wbox,8), "(<0 or 0 => contained)")
# in final cell sig is NOT tied to y by the saddle eq -- box covers sig in [0,1]; but truth has sig_true(y);
# breg_abs = max|breg| over the whole box must exceed |Btrue| for ALL true pairs AND the box corners; that's
# what the code uses (|e^{i(phi+pi/4)}-1| <= |phi+pi/4| <= (pi/rt2) y breg_abs). Check |Btrue| <= max|box|:
wabs = mp.mpf('-9')
for _ in range(400):
    y = mp.mpf(10)**random.uniform(-12, mp.log10(1.2533e-5))
    N = pi/(2*y*y); th = pi/(sqrt2*y)
    f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
    rho = mp.findroot(f, 1+y*y/2, tol=mp.mpf('1e-50'))
    sig = (rho-1)/y**2; x0=1/(rho*y); t=rho*y
    phi_true = mp.im(Phi(x0,N)) - th
    Btrue = (phi_true + pi/4)/((pi/sqrt2)*y)
    wabs = max(wabs, abs(Btrue))
print("final-cell range: max |(phi+pi/4)/((pi/rt2)y)| at true saddles =", mp.nstr(wabs,10),
      " (code's breg_abs over box [0,yJ]x[0,1]xA1box must exceed this; code prints d6 with breg_abs)")
