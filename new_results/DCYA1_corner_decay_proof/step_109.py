import sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=40
# mpmath MV0 pre (correct, from MV0 function)
a=2*mp.pi/3; N=400
pre_mp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)
print("mpmath |pre| =", float(abs(pre_mp)))
# arb pre
from flint import arb, acb, ctx
ctx.dps=60
rt2=arb(2).sqrt(); pi=arb.pi(); alpha=2*pi/3
gplus=(2*alpha.sin()).sqrt()*(acb(0,1)*(pi/4-alpha/2)).exp()
pre_arb=gplus*(acb(0,1)/pi)*(arb("1.5").gamma())*(arb("400")-arb("1.5")).gamma()/arb("400").gamma()
print("arb |pre| =", float(abs(pre_arb)))
print("ratio:", float(abs(pre_mp))/float(abs(pre_arb)))
# check gamma ratio
print("mpmath Gamma(N-1.5)/Gamma(N) =", float(mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)))
print("arb Gamma(398.5)/Gamma(400) =", float((arb("400")-arb("1.5")).gamma()/arb("400").gamma()))
