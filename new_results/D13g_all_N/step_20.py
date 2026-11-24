
from flint import arb
def iv(lo, hi):
    lo = arb(lo).lower(); hi = arb(hi).upper()
    return arb((lo+hi)/2, (hi-lo)/2)
e = iv(arb("1.0826"), arb("1.3274"))
print("iv test:", e)
mid = (arb("1.0826")+arb("1.3274"))/2
rad = (arb("1.3274")-arb("1.0826"))/2
print("mid:", mid, "rad:", rad)
print("arb(mid, rad):", arb(mid, rad))
print("arb(float(mid), float(rad)):", arb(float(mid), float(rad)))
