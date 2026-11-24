from flint import arb, acb, ctx
ctx.dps = 50
x = arb('0.5')
y = acb.log(x)   # unbound method call with arb arg — suspected crasher
print(y)
