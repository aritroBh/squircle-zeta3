import math
from flint import arb, acb, ctx
ctx.dps=60
rt2=arb(2).sqrt(); pi=arb.pi(); alpha=2*pi/3
# test sc computation
N=arb("400")
sc=N**arb("-1.5")*(-(2*alpha*N).sqrt()).exp()
print("sc:", sc)
th=(2*alpha*N).sqrt()
print("th:", th, "e^{-th}:", (-th).exp())
# test gamma
print("Gamma(398.5):", (N-arb("1.5")).gamma())
print("Gamma(1.5):", arb("1.5").gamma())
# test the model moment with interval
b=(acb(0,1)*pi/4).exp()
def gmoment(n,p2,om):
    n=int(n); P=[acb(1)]
    for kk in range(1,n+1): P.append((om*P[kk-1]+(kk-1)*P[kk-2] if kk>=2 else om*P[kk-1])/p2)
    return (-acb(0,1))**n*(2*pi/(-p2)).sqrt()*(om**2/(2*p2)).exp()*P[n]
p2=acb(-0.05,0.01); om=arb("-0.08")
print("gmoment(0):", gmoment(0,p2,om))
