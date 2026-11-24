from flint import arb, acb, ctx
ctx.dps = 50
x = arb('0.5')
z = acb.log(x)
print(z, type(z))
rr = z.real
print(rr)
