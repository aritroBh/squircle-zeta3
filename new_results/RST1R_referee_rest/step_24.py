from flint import arb, acb, ctx
ctx.dps = 60
PI = arb.pi(); SQ2 = arb(2).sqrt()
W = acb(1,1)/SQ2
r = arb(2); q = arb(3)
s = arb(1)/2 + r*W
rho = abs(s)
print("rho:", rho, type(rho))
x = rho/(rho+q/SQ2)
print("x:", x, type(x))
# try log via arb
print("arb log:", x.log())
