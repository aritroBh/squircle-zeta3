
from flint import arb, ctx
ctx.dps = 70
print("dps=70: arb(1.205, 0.1224) =", arb(1.205, 0.1224))
print("dps=70: arb(1.0826, 0.1224) =", arb(1.0826, 0.1224))
import flint
print("ctx.prec:", ctx.prec, "ctx.dps:", ctx.dps)
print("arb(1.205,arb(0.1224)):", arb(arb(1.205), arb(0.1224)))
# maybe tuple form
print("arb((1.205,),  ...)? try str:", arb("1.205 +/- 0.1224"))
