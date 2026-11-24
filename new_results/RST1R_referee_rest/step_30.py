from flint import arb, ctx
ctx.dps = 60
PI = arb.pi(); SQ2 = arb(2).sqrt()
print("ok1")
a = PI/(2*SQ2)
print("ok2", a)
