from flint import arb, acb, ctx
ctx.dps = 60
PI = arb.pi(); SQ2 = arb(2).sqrt(); a = PI/(2*SQ2)
W = acb(1,1)/SQ2
r = arb(2); q = arb(3); j = 5
s = arb(1)/2 + r*W
rho = abs(s)
print("a ok")
t1 = -a*r
print("t1 ok", t1)
z = acb.log(rho/(rho+q/SQ2))
print("log ok", z)
