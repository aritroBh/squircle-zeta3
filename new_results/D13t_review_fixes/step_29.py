import mpmath as mp
mp.mp.dps = 100
mpr2=mp.sqrt(2); mpp=mp.pi
# (4b) bracket at tiny t with 100 dps
wlo = mp.mpf('9')
for k in range(2001):
    t = mp.mpf(10)**(mp.mpf(-12) + mp.mpf(8)*k/2000)   # t in [1e-12, 1e-4]
    zz = mpr2*t/(2+mpr2*t)
    A1n = (1-mp.atan(zz)/zz)/(t*t)
    L = 2/(3*(2+mpr2*t)**2) - 4*t*t/(5*(2+mpr2*t)**4)
    wlo = min(wlo, A1n-L)
print(f"(4b) 100dps, t in [1e-12,1e-4]: min(A1n - L) = {mp.nstr(wlo,10)} (>0 => the -6.4e-38 was 60-dps noise)")

# (5) arb box containment on final cell
from flint import arb, ctx
ctx.dps = 70
art2 = arb(2).sqrt()
def iv(lo,hi):
    lo=arb(lo).lower(); hi=arb(hi).upper(); return arb((lo+hi)/2,(hi-lo)/2)
yJ = (arb.pi()/(2*arb(10)**10)).sqrt()
Tval = arb(( (1+yJ**2)*yJ ).upper())
Tbox = iv(0, Tval)
rhoB = (1 + Tbox**2/(1+art2*Tbox)).sqrt()
s1B = -Tbox/((1+art2*Tbox)*(1+rhoB))
A1nB = iv(arb(2/(3*(2+art2*Tbox)**2) - 4*Tbox**2/(5*(2+art2*Tbox)**4)).lower(),
          arb(2/(3*(2+art2*Tbox)**2)).upper())
B1box = arb("0.5") + art2*s1B + s1B**2 - rhoB**2*A1nB
Btrue_box = rhoB*B1box/(2+art2*Tbox)
blo = mp.mpf(float(arb(Btrue_box.lower()))); bhi = mp.mpf(float(arb(Btrue_box.upper())))
mp.mp.dps = 60
ma=mpp/(2*mpr2); mbb=mp.e**(1j*mpp/4)
Tv = mp.mpf(float(Tval))
viol = mp.mpf('-9'); Bt=None
for k in range(1,100001):
    t = Tv*k/100000
    rho = mp.sqrt(1+t*t/(1+mpr2*t)); y=t/rho; N=mpp/(2*y*y); th=mpp/(mpr2*y); x0=1/(rho*y)
    Bt=(mp.im((-1+1j)*ma*x0 + N*mp.log(mbb*x0/(1+mbb*x0)))-th+mpp/4)/((mpp/mpr2)*y)
    viol = max(viol, max(blo-Bt, Bt-bhi))
print(f"(5) B_true(t) in arb Btrue_box over 100000 pts t in (0,Tval]: max violation = {mp.nstr(viol,10)}")
print(f"    B_true(Tval) = {mp.nstr(Bt,15)} vs 1/6 = {mp.nstr(mp.mpf(1)/6,15)}; box = [{mp.nstr(blo,15)},{mp.nstr(bhi,15)}]")
print(f"    B_abs = 0.16666814452037215... ; breg_abs = 0.75 (old) >= B_abs: {0.75 >= float(arb(Btrue_box.upper()))}")
