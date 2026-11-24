from flint import arb, acb, ctx
ctx.dps=60
def iv(lo,hi):
    lo=arb(lo).lower(); hi=arb(hi).upper(); return arb((lo+hi)/2,(hi-lo)/2)
# Gamma(N-3/2)/Gamma(N) = B(N-3/2, 3/2)/Gamma(3/2) = (2/sqrt(pi)) B(N-3/2,3/2)
print("has beta:", hasattr(arb(2),'beta'))
N=iv("400","404")
# check beta
bval=(N-arb("1.5")).beta(arb("1.5"))
print("B(N-1.5,1.5):",bval)
rtpi=arb.pi().sqrt()
ratio=(2/rtpi)*bval
print("(2/sqrt pi) B =", ratio)
# verify vs point
import mpmath as mp
mp.mp.dps=30
print("mpmath Gamma(398.5)/Gamma(400) =", float(mp.gamma(398.5)/mp.gamma(400)))
print("interval contains?", arb("0.00012558829639990938") in ratio)
