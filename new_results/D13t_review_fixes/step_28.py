import mpmath as mp
mp.mp.dps = 60
mpr2=mp.sqrt(2); mpp=mp.pi
# (4) try to break the A1n bracket over a huge t-range (dense)
wlo = mp.mpf('9'); whi = mp.mpf('9'); wat=None
for k in range(20001):
    t = mp.mpf(10)**(mp.mpf(-12) + mp.mpf(20)*k/20000)   # t in [1e-12, 1e8]
    zz = mpr2*t/(2+mpr2*t)
    A1n = (1-mp.atan(zz)/zz)/(t*t)
    L = 2/(3*(2+mpr2*t)**2) - 4*t*t/(5*(2+mpr2*t)**4)
    U = 2/(3*(2+mpr2*t)**2)
    m1 = A1n-L; m2 = U-A1n
    if min(m1,m2) < min(wlo,whi): wat=t
    wlo=min(wlo,m1); whi=min(whi,m2)
print(f"(4) A1n bracket over t in [1e-12,1e8], 20001 pts: min(A1n-L) = {mp.nstr(wlo,8)}, min(U-A1n) = {mp.nstr(whi,8)} (>=0 => never fails) at t={mp.nstr(wat,6)}")

# (5) arb box containment on the final cell (dense)
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
print(f"    A1nB = {A1nB} (well-formed: radius>0 = {arb(A1nB.rad())>0})")
B1box = arb("0.5") + art2*s1B + s1B**2 - rhoB**2*A1nB
Btrue_box = rhoB*B1box/(2+art2*Tbox)
B_abs = max(abs(arb(Btrue_box.lower())), abs(arb(Btrue_box.upper())))
print(f"    Btrue_box = {Btrue_box}, B_abs = {B_abs}")
blo, bhi = mp.mpf(str(arb(Btrue_box.lower()))), mp.mpf(str(arb(Btrue_box.upper())))
viol = mp.mpf('-9'); Blimit=None
for k in range(1,100001):
    t = mp.mpf(str(Tval))*k/100000
    rho = mp.sqrt(1+t*t/(1+mpr2*t)); y=t/rho; N=mpp/(2*y*y); th=mpp/(mpr2*y); x0=1/(rho*y)
    import mpmath as m2
    ma=mpp/(2*mpr2); mbb=m2.e**(1j*mpp/4)
    Bt=(m2.im((-1+1j)*ma*x0 + N*m2.log(mbb*x0/(1+mbb*x0)))-th+mpp/4)/((mpp/mpr2)*y)
    viol = max(viol, max(blo-Bt, Bt-bhi))
print(f"    B_true(t) in arb box over 100000 pts t in (0,Tval]: max violation = {mp.nstr(viol,8)} (<=0 contained)")
print(f"    B_true -> 1/6 as t->0? B_true(Tval) = {mp.nstr(Bt,15)}, 1/6 = {mp.nstr(mpp.mpf(1)/6,15)}")
