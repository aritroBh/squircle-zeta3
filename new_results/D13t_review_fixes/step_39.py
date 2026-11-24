import sys, math
sys.path.insert(0, 'research/08_theory/D13s_fix_recertify/D13s_work/fixed')
from flint import arb, acb, ctx
ctx.dps = 70
import importlib
st = importlib.import_module('step_38')   # fixed version (runs its __main__? no, guarded)
iv, pi, rt2, h3encl = st.iv, st.pi, st.rt2, st.h3encl
import mpmath as mp
mp.mp.dps = 80
r2=mp.sqrt(2); mpp=mp.pi; ma=mpp/(2*r2); mbb=mp.e**(1j*mpp/4)
def lhs_true(yv, sgv):
    rho = 1+yv*yv*sgv; t=rho*yv; N=mpp/(2*yv*yv); th=mpp/(r2*yv); x0=1/(rho*yv)
    return mp.re((-1+1j)*ma*x0 + N*mp.log(mbb*x0/(1+mbb*x0))) + th

def mS_box(y1,y2,s1,s2):
    y = iv(y1,y2); sig = iv(s1,s2)
    rho = 1+y*y*sig; y0 = rho*y; N = pi/(2*y*y)
    return (N*y0**3/2)*h3encl(y0)

# first certified cells (ratio 1.05) : [100,105], [105,110.25]
for (n1,n2) in [(100.,105.),(105.,110.25),(120.,126.),(130.,136.5)]:
    y1 = math.sqrt(math.pi/(2*n2)); y2 = math.sqrt(math.pi/(2*n1))
    # subcell sigma in [0.75,1]
    mS = mS_box(y1,y2,0.75,1.0)
    mlo = mp.mpf(float(arb(mS.lower())))
    # true lhs at (y2, sigma=1) (N=n1 corner)
    lt = lhs_true(mp.mpf(y2), mp.mpf(1))
    print(f"cell N in [{n1},{n2}], sig [0.75,1]: mS.lower = {mp.nstr(mlo,10)}, lhs(y2,sig=1) = {mp.nstr(lt,10)}, CONTAINED? {lt >= mlo}  (margin {mp.nstr(lt-mlo,6)})")
# and at the saddle for the same cell
for (n1,n2) in [(100.,105.)]:
    y1 = math.sqrt(math.pi/(2*n2)); y2 = math.sqrt(math.pi/(2*n1))
    mS = mS_box(y1,y2,0.25,0.5)
    mlo = mp.mpf(float(arb(mS.lower()))); mhi = mp.mpf(float(arb(mS.upper())))
    yv = mp.sqrt(mpp/(2*mp.mpf(n1)))
    f = lambda r: r**2*(1+r2*r*yv)/(1+r2*r*yv+r**2*yv**2) - 1
    rho = mp.findroot(f, 1+yv*yv/2); sg = (rho-1)/yv**2
    lt = lhs_true(yv, sg)
    print(f"cell N in [{n1},{n2}], saddle subcell: sigma* = {mp.nstr(sg,8)}, lhs = {mp.nstr(lt,10)} in [{mp.nstr(mlo,10)},{mp.nstr(mhi,10)}]? {mlo <= lt <= mhi} (margins {mp.nstr(lt-mlo,6)}, {mp.nstr(mhi-lt,6)})")
