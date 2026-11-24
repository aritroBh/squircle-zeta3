from flint import arb, acb, ctx
ctx.dps=60
def iv(lo,hi):
    lo=arb(lo).lower(); hi=arb(hi).upper(); return arb((lo+hi)/2,(hi-lo)/2)
N=iv("400","404")
# check lgamma availability
print("has lgamma:", hasattr(arb(2),'lgamma'))
try:
    print("(N-1.5).lgamma():",(N-arb("1.5")).lgamma())
except Exception as e:
    print("err:",e)
try:
    print("N.lgamma():",N.lgamma())
except Exception as e:
    print("err:",e)
# compute ratio via lgamma
lg1=(N-arb("1.5")).lgamma(); lg2=N.lgamma()
print("lgamma ratio exp:", (lg1-lg2).exp())
