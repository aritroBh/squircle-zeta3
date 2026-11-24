import sys, math
sys.path.insert(0, '<project>/research/08_theory/D13s_fix_recertify/D13s_work/fixed')
from flint import arb, ctx
ctx.dps = 70
import importlib
st = importlib.import_module('step_38')
iv, h3encl, pi, rt2 = st.iv, st.h3encl, st.pi, st.rt2
sq2 = rt2; a2 = pi*sq2/4
def truth_pt(y, sig):
    y = arb(y); sig = arb(sig)
    rho = 1 + y*y*sig; t = rho*y; N = pi/(2*y*y)
    w = 1 - (1+sq2*t+t*t)*(-sq2*t).exp(); psi = -(1-w).log()
    return (N/2)*psi - (a2/t)*(rho-1)**2
def coded_mS(y1,y2,s1,s2):
    y = iv(y1,y2); sig = iv(s1,s2)
    rho = 1 + y*y*sig; y0 = rho*y; N = pi/(2*y*y)
    return (N*y0**3/2)*h3encl(y0)

print("cell-by-cell at corner (y=y2 i.e. N=n1, sigma=1): TIGHT box vs CODED box")
print(f"{'cell N':>18} {'truth':>10} {'tight_lo':>10} {'tight marg':>11} {'coded_lo':>10} {'coded marg':>11} {'slack factors':>22}")
for (n1,n2) in [(100.,105.),(105.,110.25),(120.,126.),(130.,136.5),(373.3456322341576,392.0129138458655),(392.0129138458655,400.)]:
    y1 = math.sqrt(math.pi/(2*n2)); y2 = math.sqrt(math.pi/(2*n1))
    yv = arb(y2); t = yv*(1+yv*yv); N = pi/(2*yv*yv)
    tight_lo = (N*t**3/2)*(sq2/3)/(1+sq2*t+t*t)
    tv = truth_pt(y2, 1.0)
    mS = coded_mS(y1,y2,0.75,1.0); clo = arb(mS.lower())
    # slack decomposition: coded_lo / tight_lo
    yb = iv(y1,y2); sb = iv(0.75,1.0); rhob = 1+yb*yb*sb; y0b = rhob*yb; Nb = pi/(2*yb*yb); T = arb(y0b.upper())
    f1 = arb(Nb.lower())/N            # N.lo/N
    f2 = (arb(y0b.lower())/t)**3      # (y0.lo/t)^3
    f3 = (1+sq2*t+t*t)/(1+sq2*T+T*T)  # h3lo denominator T-slack
    print(f"[{n1:8.3f},{n2:8.3f}] {float(arb(tv.mid())):10.6f} {float(arb(tight_lo.mid())):10.6f} {float(arb(arb(tv.lower())-arb(tight_lo.upper()))):+11.6f} {float(clo):10.6f} {float(arb(arb(tv.lower())-clo)):+11.6f}  {float(f1):.4f}x{float(f2):.4f}x{float(f3):.4f}")
print()
print("at the saddle (cell [100,105], y=y2, sigma*=0.42857958...):")
import mpmath as mp
mp.mp.dps=50; r2=mp.sqrt(2)
y2 = math.sqrt(math.pi/(2*100.))
f=lambda tt: tt/mp.sqrt(1+tt**2/(1+r2*tt))-y2
tm=mp.findroot(f, y2*(1+y2**2/2)); rhom=mp.sqrt(1+tm**2/(1+r2*tm)); sgm=float((rhom-1)/y2**2)
yv=arb(y2); tv = truth_pt(y2, sgm); t = arb(tm); N = pi/(2*yv*yv)
tight_lo = (N*t**3/2)*(sq2/3)/(1+sq2*t+t*t)
mS = coded_mS(math.sqrt(math.pi/(2*105.)), y2, 0.25, 0.5)
print(f"  truth(saddle) = {float(arb(tv.mid())):.10f}")
print(f"  tight_lo = {float(arb(tight_lo.mid())):.10f}, tight margin = {float(arb(arb(tv.lower())-arb(tight_lo.upper()))):+.6e}")
print(f"  coded mS = [{float(arb(mS.lower())):.10f}, {float(arb(mS.upper())):.10f}], coded margins = +{float(arb(arb(tv.lower())-arb(mS.lower()))):.6f} / +{float(arb(arb(mS.upper())-arb(tv.upper()))):.6f}")
