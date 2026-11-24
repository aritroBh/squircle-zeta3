import math
from flint import arb, acb, ctx
ctx.dps=60
rt2=arb(2).sqrt(); pi=arb.pi(); b=(acb(0,1)*pi/4).exp(); ialpha=acb(0,2*pi/3); alpha=2*pi/3
def iv(lo,hi):
    lo=arb(lo).lower(); hi=arb(hi).upper(); return arb((lo+hi)/2,(hi-lo)/2)
gplus=(2*alpha.sin()).sqrt()*(acb(0,1)*(pi/4-alpha/2)).exp()
def gamma_ratio_interval(Ncell):
    # Gamma(N-3/2)/Gamma(N) via lgamma at center + digamma variation
    n1=arb(Ncell).lower(); n2=arb(Ncell).upper(); nc=(n1+n2)/2
    lg1=(nc-arb("1.5")).lgamma(); lg2=nc.lgamma()
    logr_c=lg1-lg2
    # |d log r/dN| <= digamma(N-1.5)+|digamma(N)| ~ 2 log N bound; use 1.6/N1 conservatively... use actual bound
    # digamma(N-1.5)-digamma(N) = sum_{j>=0}[1/(N+j)-1/(N-1.5+j)], |.|<=1.5/(N-1.5)/1 ... use 1.6/(n1-1.5)
    dn1=arb(n1).lower()
    dg=arb("1.6")/(dn1-arb("1.5"))
    dlogr=dg*(n2-n1)
    return iv(logr_c-dlogr, logr_c+dlogr).exp()
Ncell=iv("400","404")
gr=gamma_ratio_interval(Ncell)
print("gamma ratio interval:",gr)
print("contains point 0.00012558829...?", arb("0.00012558829639990938") in gr)
# now prefactor
def preF(Ncell):
    return gplus*(acb(0,1)/pi)*(arb("1.5").gamma())*gamma_ratio_interval(Ncell)
pre=preF(Ncell)
print("pre:",pre)
print("abs(pre):",abs(pre))
