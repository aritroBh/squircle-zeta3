from flint import arb, acb, ctx
ctx.dps = 30
u = arb("0.6","0.7")
print("arb('0.6','0.7') =", u)  # check meaning
u2 = arb(0.6, 0.7)
print("arb(0.6,0.7) =", u2)
u3 = arb(0.6)
print("arb(0.6) =", u3)
# interval [lo,hi]: use arb(lo, hi)?  try arb(lo, rad=(hi-lo)/2)
lo = arb("0.6"); hi = arb("0.7")
ball = lo + (hi-lo)/2 * arb(-1,1)
print("proper interval:", ball)
print("acb of it:", acb(ball))
