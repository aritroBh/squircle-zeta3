import os
os.chdir('<project>/research/08_theory/RST1_rest_exponent')
step = r'''"""RST1 step_12: detailed q2 constant decomposition at k=99 and monotonicity check."""
import mpmath as mp
mp.mp.dps=50
from mpmath import mpf, pi, sqrt, exp
a=pi/(2*sqrt(2)); sq2=sqrt(2); sqpi=sqrt(pi); sq2pi=sqrt(2*pi)
CA=(4*sq2)**mpf('1.5')*sqrt(125/(2*sq2)); CL=2**mpf('1.5')/(3*a); CR=1/a
print("CA=",mp.nstr(CA,10)," CL=",mp.nstr(CL,10)," CR=",mp.nstr(CR,10)," a=",mp.nstr(a,10))
def FHg(j,q):
    r_=sqrt(2*j*q/pi); c0=1+q/sq2; spj=sqrt(pi*j*q)
    return sqpi/q*(CL*exp(-spj/4)+CA/r_+CR*(sq2*r_)**mpf('-1.5'))*exp(a*c0/2)*exp(-spj)
def tg(k):
    N=k+1; return N**mpf('-0.5')*exp(-sq2pi*sqrt(N))
k=99
i0=FHg(k,2); i1=2*FHg(k-1,2)
q2=(8/pi)*(i0+i1)
print("q2 bound at k=99:",mp.nstr(q2,10))
print("q2/target at k=99:",mp.nstr(q2/tg(k),10))
# sub-pieces normalized for i=0 and i=1
for nm,(j,q) in [("i0",(k,2)),("i1",(k-1,2))]:
    r_=sqrt(2*j*q/pi); c0=1+q/sq2; spj=sqrt(pi*j*q)
    fac=(8/pi)*(2 if nm=="i1" else 1)*sqpi/q*exp(a*c0/2)*exp(-spj)/tg(k)
    le=fac*CL*exp(-spj/4); mi=fac*CA/r_; ri=fac*CR*(sq2*r_)**mpf('-1.5')
    print(f"{nm}: left={mp.nstr(le,8)} mid={mp.nstr(mi,8)} right={mp.nstr(ri,8)}  (normalized to target)")
# monotonicity of each normalized sub-piece for k>=99: check derivative sign numerically over grid
def subpieces(k):
    out=[]
    for w,(j,q) in [(1,(k,2)),(2,(k-1,2))]:
        r_=sqrt(2*j*q/pi); c0=1+q/sq2; spj=sqrt(pi*j*q)
        fac=(8/pi)*w*sqpi/q*exp(a*c0/2)*exp(-spj)/tg(k)
        out += [fac*CL*exp(-spj/4), fac*CA/r_, fac*CR*(sq2*r_)**mpf('-1.5')]
    return out
import numpy as np
bad=0
for idx in range(6):
    for k in range(99,2000):
        if subpieces(k+1)[idx] > subpieces(k)[idx]*(1+1e-12):
            bad+=1
            if bad<5: print("INCREASE at piece",idx,"k",k)
print("monotonicity violations (k in [99,2000)):",bad)
# sum of all six at k=99:
print("sum of six q2 sub-pieces at k=99:", mp.nstr(sum(subpieces(99)),10))
'''
open('step_12.py','w').write(step)
import subprocess
r = subprocess.run(['bash','-c','cd <project>/research/08_theory/RST1_rest_exponent && nice <project>/.venv/bin/python step_12.py'], capture_output=True, text=True, timeout=2400)
print(r.stdout)
print(r.stderr[-1500:] if r.stderr else '')
