from flint import arb, acb, ctx
ctx.dps = 60
PI = arb.pi(); SQ2 = arb(2).sqrt(); a = PI/(2*SQ2)
W = acb(1,1)/SQ2
r = arb(2); q = arb(3); j = 5
s = arb(1)/2 + r*W
rho = abs(s)
print(type(rho))
exact = -a*r + j*acb.log(rho/(rho+q/SQ2)).real
print(exact)
Etil = -a*r - SQ2*j*q/(2*r+1+q/SQ2)
print(Etil - exact)
