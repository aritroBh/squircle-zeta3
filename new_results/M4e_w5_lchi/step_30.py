import time, math, json
import mpmath as mp
# check flint acb Hurwitz zeta availability
from flint import acb, arb, ctx
print("acb.zeta signature test:")
ctx.prec = 300
try:
    z = acb.zeta(acb(5), acb(1)/8)
    print("two-arg works:", z)
except Exception as e:
    print("two-arg failed:", e)
