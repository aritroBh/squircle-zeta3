from flint import arb, acb, ctx
ctx.dps=30
z=acb(1,2)
print("real:", z.real, type(z.real))
print("imag:", z.imag, type(z.imag))
# test .real() vs .real
