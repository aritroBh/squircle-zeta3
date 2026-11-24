import math
from flint import arb, acb, ctx
ctx.dps=60
rt2=arb(2).sqrt(); pi=arb.pi(); alpha=2*pi/3
def iv(lo,hi):
    lo=arb(lo).lower(); hi=arb(hi).upper(); return arb((lo+hi)/2,(hi-lo)/2)
# old (amplifying) RePhi
def RePhi_old(t,N): return -t/rt2+(N-1)*t.log()-(N-arb("1.5"))/2*(t**2+rt2*alpha*t+alpha**2).log()
# new (combined) RePhi
def RePhi(t,N):
    u=rt2*alpha/t + alpha**2/t**2
    return -t/rt2 + arb("0.5")*t.log() - (N-arb("1.5"))/2*(1+u).log()
# verify match and tightness
import mpmath as mp
mp.mp.dps=30
bmp=mp.e**(1j*mp.pi/4)
def RePhi_m(t,N): return -t/mp.sqrt(2)+(N-1)*mp.log(t)-(N-mp.mpf('1.5'))/2*mp.log(t**2+mp.sqrt(2)*2.0944*t+2.0944**2)
for N in [400,1600]:
    t=(arb(N)*alpha).sqrt()
    print(f"N={N}: combined={float(arb(RePhi(t,N)).mid()):.6f} mpmath={float(RePhi_m(float(t),N)):.6f}")
# check interval tightness
Ncell=iv("400","404")
# t0 range ~ [29.178,29.896]
t=iv("29.178","29.896")
print("RePhi interval on t-range:", RePhi(t,Ncell))
print("old RePhi interval:", RePhi_old(t,Ncell))
