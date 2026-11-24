
import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); pi = mp.pi; a = pi/(2*sqrt2); b = mp.e**(1j*pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(b*x/(1+b*x))

# Exact identity: RePhi(x0)+th = -(N/2)S(t) - (pi sqrt2/(4 t))(rho-1)^2, t = rho y, saddle branch rho^2 = 1 + t^2/(1+sqrt2 t)
worst = mp.mpf('0')
for k in range(-400, 401):
    t = mp.mpf(10)**(mp.mpf(k)/40)   # t in [1e-10, 1e10]
    rho = mp.sqrt(1 + t*t/(1+sqrt2*t))
    y = t/rho; N = pi/(2*y*y); th = pi/(sqrt2*y); x0 = 1/(rho*y)
    S = mp.log(1+sqrt2*t+t*t) - sqrt2*t
    lhs = mp.re(Phi(x0,N)) + th
    rhs = -(N/2)*S - (pi*sqrt2/(4*t))*(rho-1)**2
    worst = max(worst, abs(lhs-rhs)/abs(lhs))
print("H3 identity RePhi0+th = -(N/2)S - (pi rt2/(4t))(rho-1)^2 : worst rel resid =", mp.nstr(worst, 8))

# h3-box containment: box = [(N t^3/2)(rt2/3)/(1+rt2 t+t^2), (N t^3/2)(rt2/3)]; truth inside for all t>0 (saddle branch)
minlo = mp.mpf('9'); minhi = mp.mpf('9')
for k in range(-400, 401):
    t = mp.mpf(10)**(mp.mpf(k)/40)
    rho = mp.sqrt(1 + t*t/(1+sqrt2*t))
    y = t/rho; N = pi/(2*y*y); th = pi/(sqrt2*y); x0 = 1/(rho*y)
    truth = mp.re(Phi(x0,N)) + th
    blo = (N*t**3/2)*(sqrt2/3)/(1+sqrt2*t+t*t)
    bhi = (N*t**3/2)*(sqrt2/3)
    minlo = min(minlo, truth - blo)
    minhi = min(minhi, bhi - truth)
print("h3-box containment: min(truth - box_lo) =", mp.nstr(minlo,8), "  min(box_hi - truth) =", mp.nstr(minhi,8))

# psi lemma check: psi(t)/t^3 = (rt2 t - log(1+rt2 t+t^2))/t^3 in [rt2/3/(1+rt2 t+t^2), rt2/3]
mn = mp.mpf('9'); mx = mp.mpf('-9'); mn2 = mp.mpf('9')
for k in range(-400, 401):
    t = mp.mpf(10)**(mp.mpf(k)/40)
    psi = sqrt2*t - mp.log(1+sqrt2*t+t*t)
    r = psi/t**3
    mn = min(mn, r - (sqrt2/3)/(1+sqrt2*t+t*t))
    mx = max(mx, r - sqrt2/3)
print("psi/t^3 - rt2/3/(1+..): min =", mp.nstr(mn,8), " ; psi/t^3 - rt2/3: max =", mp.nstr(mx,8))

# H1: exact Btrue regular form  Btrue = breg_h[A1'] - 1/12
worst2 = mp.mpf('0'); mn12 = mp.mpf('9')
for _ in range(400):
    import random; random.seed(_)
    t = mp.mpf(10)**random.uniform(-9, 1)
    rho = mp.sqrt(1 + t*t/(1+sqrt2*t))
    y = t/rho; N = pi/(2*y*y); th = pi/(sqrt2*y); x0 = 1/(rho*y)
    sig = (rho-1)/y**2
    phi_true = mp.im(Phi(x0,N)) - th
    zz = sqrt2*t/(2+sqrt2*t)
    Ap = 1 - mp.atan(zz)/zz
    breg_h = -sig/(2*rho) + (-2*y*sig + 2*sqrt2*sig + sqrt2*rho - 2*sqrt2*rho*Ap)/(2*sqrt2*(2+sqrt2*t))
    Btrue = (phi_true + pi/4)/((pi/sqrt2)*y)
    worst2 = max(worst2, abs(Btrue - (breg_h - mp.mpf(1)/12)))
    mn12 = min(mn12, breg_h - Btrue)
print("H1: Btrue = breg_h - 1/12: worst resid =", mp.nstr(worst2,8), " ; min(breg_h - Btrue) =", mp.nstr(mn12, 12))
