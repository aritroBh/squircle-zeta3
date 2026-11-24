import time
from flint import arb, acb, ctx
ctx.dps=50
# test arb basics
x=arb("2.5"); y=acb(1,2)
print("arb:", x, "acb:", y, "gamma:", arb("2.5").gamma())
print("exp:", x.exp(), "sqrt:", x.sqrt(), "log:", x.log())
print("pi:", arb.pi())
# test interval
def iv(lo,hi):
    lo=arb(lo).lower(); hi=arb(hi).upper(); return arb((lo+hi)/2,(hi-lo)/2)
z=iv("1.0","1.2")
print("iv:", z, "contains 1.1?", arb("1.1") in z)
print("OK, flint works")
