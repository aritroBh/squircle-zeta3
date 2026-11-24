import math
from flint import arb, acb, ctx
ctx.dps=60
rt2=arb(2).sqrt(); pi=arb.pi(); alpha=2*pi/3
def RePhip2(t,N):
    # cancellation-free: RePhip = -1/rt2 + num/(t q), num = (1/2)t^2 + alpha((N-1)/rt2+1/(2 rt2))t + (N-1)alpha^2
    q=t**2+rt2*alpha*t+alpha**2
    num=arb("0.5")*t**2 + alpha*((N-1)/rt2 + arb("0.5")/rt2)*t + (N-1)*alpha**2
    return -1/rt2 + num/(t*q)
# check vs direct
def RePhip(t,N):
    q=t**2+rt2*alpha*t+alpha**2
    return -1/rt2+(N-1)/t-(N-arb("1.5"))*(t+alpha/rt2)/q
for N in [400,1600]:
    for mult in [1.1,1.3,1.5,2.0]:
        t=mult*(arb(N)*alpha).sqrt()
        print(f"N={N} t={mult}sqrt: direct={float(arb(RePhip(t,N)).mid()):.4f} combined={float(arb(RePhip2(t,N)).mid()):.4f}")
