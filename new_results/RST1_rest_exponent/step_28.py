import os
os.chdir('<project>/research/08_theory/RST1_rest_exponent')
step = r'''"""RST1 step_13: exact-form F_H lemma -> smaller C_q2.
mid piece uses x_b=(r_*-c0)/(2 sqrt2) and c2=8 sqrt2 j q/(4 r_*+c0)^3 exactly (no relaxation)."""
import mpmath as mp
mp.mp.dps=50
from mpmath import mpf, pi, sqrt, exp
a=pi/(2*sqrt(2)); sq2=sqrt(2); sqpi=sqrt(pi); sq2pi=sqrt(2*pi)
CL=2**mpf('1.5')/(3*a); CR=1/a
def FHg_exact(j,q):
    r_=sqrt(2*j*q/pi); c0=1+q/sq2; spj=sqrt(pi*j*q)
    rb=(r_-c0)/2; xb=mpf('0.5')+rb/sq2
    c2=8*sq2*j*q/(4*r_+c0)**3
    mid=xb**mpf('-1.5')*sqrt(2*pi/c2)
    return sqpi/q*(CL*exp(-spj/4)+mid+CR*(sq2*r_)**mpf('-1.5'))*exp(a*c0/2)*exp(-spj)
def FHg_relaxed(j,q):
    r_=sqrt(2*j*q/pi); c0=1+q/sq2; spj=sqrt(pi*j*q)
    CA=(4*sq2)**mpf('1.5')*sqrt(125/(2*sq2))
    return sqpi/q*(CL*exp(-spj/4)+CA/r_+CR*(sq2*r_)**mpf('-1.5'))*exp(a*c0/2)*exp(-spj)
def tg(k):
    N=k+1; return N**mpf('-0.5')*exp(-sq2pi*sqrt(N))
for k in [99,400,10000]:
    q2e=(8/pi)*(FHg_exact(k,2)+2*FHg_exact(k-1,2))/tg(k)
    q2r=(8/pi)*(FHg_relaxed(k,2)+2*FHg_relaxed(k-1,2))/tg(k)
    print(f"k={k}: exact C_q2={mp.nstr(q2e,8)}  relaxed C_q2={mp.nstr(q2r,8)}")
# monotonicity of exact mid sub-pieces on grid
def mids_exact(k):
    out=[]
    for w,(j,q) in [(1,(k,2)),(2,(k-1,2))]:
        r_=sqrt(2*j*q/pi); c0=1+q/sq2; spj=sqrt(pi*j*q)
        rb=(r_-c0)/2; xb=mpf('0.5')+rb/sq2
        c2=8*sq2*j*q/(4*r_+c0)**3
        mid=xb**mpf('-1.5')*sqrt(2*pi/c2)
        fac=(8/pi)*w*sqpi/q*exp(a*c0/2)*exp(-spj)/tg(k)
        out += [fac*CL*exp(-spj/4), fac*mid, fac*CR*(sq2*r_)**mpf('-1.5')]
    return out
bad=0
for idx in range(6):
    for k in range(99,3000):
        if mids_exact(k+1)[idx] > mids_exact(k)[idx]*(1+1e-12): bad+=1
print("exact-form monotonicity violations on [99,3000):",bad)
print("exact sum at 99:", mp.nstr(sum(mids_exact(99)),8))
'''
open('step_13.py','w').write(step)
import subprocess
r = subprocess.run(['bash','-c','cd <project>/research/08_theory/RST1_rest_exponent && nice <project>/.venv/bin/python step_13.py'], capture_output=True, text=True, timeout=2400)
print(r.stdout)
print(r.stderr[-1500:] if r.stderr else '')
