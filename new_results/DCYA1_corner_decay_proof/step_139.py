import math
from flint import arb, acb, ctx
ctx.dps=60
rt2=arb(2).sqrt(); pi=arb.pi(); alpha=2*pi/3
def iv(lo,hi):
    lo=arb(lo).lower(); hi=arb(hi).upper(); return arb((lo+hi)/2,(hi-lo)/2)
# combined c2 (all-positive numerator)
def c2_combined(t,N):
    q=t**2+rt2*alpha*t+alpha**2
    num=arb("0.5")*t**4 + rt2*alpha*(N-arb("0.5"))*t**3 + 4*alpha**2*(N-1)*t**2 + 2*rt2*alpha**3*(N-1)*t + alpha**4*(N-1)
    return num/(t**2*q**2)
def RePhipp(t,N):
    q=t**2+rt2*alpha*t+alpha**2
    return -(N-1)/t**2+(N-arb("1.5"))*(t**2+rt2*alpha*t)/q**2
def c2_old(t,N): return -RePhipp(t,N)
# verify match and positivity
import mpmath as mp
mp.mp.dps=30
def c2_m(t,N):
    q=t**2+mp.sqrt(2)*2.0944*t+2.0944**2
    return -(-(N-1)/t**2+(N-mp.mpf('1.5'))*(t**2+mp.sqrt(2)*2.0944*t)/q**2)
for N in [400,1600]:
    t=(arb(N)*alpha).sqrt()
    print(f"N={N}: combined={float(arb(c2_combined(t,N)).mid()):.6f} mpmath={float(c2_m(float(t),N)):.6f}")
# positivity on interval
for (n1,n2) in [(400,404),(1600,1616)]:
    Ncell=iv(str(n1),str(n2))
    t=iv(str((arb(n1)*alpha).sqrt().lower()),str((arb(n2)*alpha).sqrt().upper()*1.01))
    c=c2_combined(t,Ncell)
    print(f"[{n1},{n2}]: c2=[{float(arb(c).lower()):.5f},{float(arb(c).upper()):.5f}] positive? {arb(c).lower()>0}")
