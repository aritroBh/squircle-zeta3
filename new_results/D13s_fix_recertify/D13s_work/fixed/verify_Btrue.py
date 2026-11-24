
# H1: EXACT regular form  Btrue = rho*B1/(2+rt2 t)  vs direct  Btrue = (Im Phi(x0)-th+pi/4)/((pi/rt2)y)
import mpmath as mp, random
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); pi = mp.pi; a = pi/(2*sqrt2); b = mp.e**(1j*pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(b*x/(1+b*x))
def Bform(t):
    rho = mp.sqrt(1 + t*t/(1+sqrt2*t))
    s1 = -t/((1+sqrt2*t)*(1+rho))
    zz = sqrt2*t/(2+sqrt2*t)
    Ap = 1 - mp.atan(zz)/zz            # exact A1'
    A1n = Ap/(t*t)                     # exact A1n = A1'/t^2
    B1 = mp.mpf('0.5') + sqrt2*s1 + s1*s1 - rho*rho*A1n
    return rho*B1/(2+sqrt2*t)
worst = mp.mpf('0'); random.seed(7)
for _ in range(600):
    t = mp.mpf(10)**random.uniform(-12, 0)
    rho = mp.sqrt(1 + t*t/(1+sqrt2*t)); y = t/rho; N = pi/(2*y*y); th = pi/(sqrt2*y); x0 = 1/(rho*y)
    Bt = (mp.im(Phi(x0,N)) - th + pi/4)/((pi/sqrt2)*y)
    worst = max(worst, abs(Bt - Bform(t))/abs(Bt))
print("H1 exact form (exact A1n): worst rel resid over t in [1e-12,1] (600 pts) =", mp.nstr(worst, 8))

# containment in the arb Btrue_box used by the fixed final cell
from flint import arb, ctx
ctx.dps = 60
art2 = arb(2).sqrt()
def iv(lo, hi):
    lo = arb(lo).lower(); hi = arb(hi).upper()
    return arb((lo+hi)/2, (hi-lo)/2)
yJ = (arb.pi()/(2*arb(10)**10)).sqrt()
Tval = arb(( (1+yJ**2)*yJ ).upper())
Tbox = iv(0, Tval)
rhoB = (1 + Tbox**2/(1+art2*Tbox)).sqrt()
s1B = -Tbox/((1+art2*Tbox)*(1+rhoB))
A1nB = iv(arb(2/(3*(2+art2*Tbox)**2) - 4*Tbox**2/(5*(2+art2*Tbox)**4)).lower(),
          arb(2/(3*(2+art2*Tbox)**2)).upper())
B1box = arb("0.5") + art2*s1B + s1B**2 - rhoB**2*A1nB
Btrue_box = rhoB*B1box/(2+art2*Tbox)
print("arb Btrue_box =", Btrue_box)
blo, bhi = float(arb(Btrue_box.lower())), float(arb(Btrue_box.upper()))
viol = mp.mpf('-9')
for k in range(20001):
    t = float(Tval)*k/20000
    if t == 0: continue
    v = Bform(mp.mpf(t))
    viol = max(viol, max(blo - v, v - bhi))
print("Btrue(t) inside arb box over 20000 grid pts: max violation =", mp.nstr(viol, 8), "(<=0 => contained)")
