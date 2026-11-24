import time, sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=80
def gamma_m(c, M):
    Q=[mp.mpf(0)]*(M+1); Q[0]=1
    if M>=1: Q[1]=-2*c
    if M>=2: Q[2]=1
    b=[mp.sqrt(Q[0])]
    for n in range(1,M+1):
        s=Q[n]-mp.fsum(b[k]*b[n-k] for k in range(1,n))
        b.append(s/(2*b[0]))
    return b
def beta_j(u, J):
    out=[mp.mpf(1)]
    for j in range(1,J+1): out.append(-(u+j)/j*out[-1])
    return out
c=mp.cos(mp.pi/2); u=mp.mpf('-0.9'); g0=mp.mpf(1)/64
gm=gamma_m(c,30); bj=beta_j(u,30)
ana=mp.mpf(0)
for j in range(31):
    for m in range(31):
        if (j,m)==(0,0): continue
        nu=j+(m-1)*(1-u); ana+=bj[j]*gm[m]*g0**(nu+1)/(nu+1)
def F(g):
    return ((1+g)**(-1-u)*mp.sqrt(1+g**(2-2*u)-2*c*g**(1-u))-1)*g**(u-1)
# numerical part via mp.quad (default)
Iq=mp.quad(F, [g0, mp.mpf(1)/8, mp.mpf(1)/2, 1])
# manual fine: many subintervals
import numpy as np
pts=[g0]
for i in range(1,40):
    pts.append(g0+(1-g0)*i/39)
Ifine=mp.quad(F, pts)
print("ana=",mp.nstr(ana,30))
print("mp.quad num   =",mp.nstr(Iq,30))
print("manual 40-num =",mp.nstr(Ifine,30))
print("diff in num part:",mp.nstr(abs(Iq-Ifine),4))
print("M (mp.quad) =",mp.nstr(ana+Iq,32))
print("M (manual)  =",mp.nstr(ana+Ifine,32))
