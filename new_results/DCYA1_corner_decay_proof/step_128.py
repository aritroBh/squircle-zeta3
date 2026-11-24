from flint import arb, acb, ctx
ctx.dps=60
N=arb("400",10)  # wait, arb("400",10) is mid+rad; let me use iv
def iv(lo,hi):
    lo=arb(lo).lower(); hi=arb(hi).upper(); return arb((lo+hi)/2,(hi-lo)/2)
N=iv("400","404")
print("N:",N)
print("(N-1.5):",(N-arb("1.5")))
print("(N-1.5).gamma():",(N-arb("1.5")).gamma())
print("N.gamma():",N.gamma())
print("ratio:",(N-arb("1.5")).gamma()/N.gamma())
print("g15:",arb("1.5").gamma())
