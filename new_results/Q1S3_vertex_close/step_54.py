import importlib, q3_engine, q1s_exact, q3_eval
for m in (q1s_exact, q3_engine, q3_eval): importlib.reload(m)
from q3_engine import *
from q3_eval import GEval
from q1s_exact import femb
import mpmath as mp
mp.mp.dps=80
r2=mp.sqrt(2); a=1+r2; b=r2-1
ge=GEval(80, M=48, gln=5)
def acc_value(acc):
    tot=mp.mpc(0)
    for (w,lvm),c in acc.W.items():
        wn=tuple(femb(s,1) for s in w)
        v=ge.word(wn)
        for f in lvm: v*=mp.log(femb(f[1],1))
        tot+=femb(c,1)*v
    for lvm,c in acc.P.items():
        v=mp.mpc(1)
        for f in lvm: v*=mp.log(femb(f[1],1))
        tot+=femb(c,1)*v
    tot+=femb(acc.R,1)
    return tot
def hg1(y): return mp.log((1-b*y)*(1+a*y)/(2*y))
def hg2(y): return mp.log((1+a*a*y*y)*(1+b*b*y*y)/((1-b*y)**2*(1+a*y)**2))
P12 = le_prod(g1,g2)
# kernel (1,0):
acc0=Acc(); int_simple_kernel(acc0, ONE, s0, P12)
d0 = mp.quad(lambda y: hg1(y)*hg2(y)/y, [0, mp.mpf('0.4'), 1])
print("kernel 1/y: engine", mp.nstr(acc_value(acc0),18), " direct", mp.nstr(d0,18))
# kernel (sqrt2/2, a):
hr2 = fscale(r2_, Fr(1,2))
accA=Acc(); int_simple_kernel(accA, hr2, sa, P12)
dA = mp.quad(lambda y: (r2/2)/(y-a)*hg1(y)*hg2(y), [0,1])
print("kernel (r2/2)/(y-a): engine", mp.nstr(acc_value(accA),18), " direct", mp.nstr(dA,18))
# kernel (-sqrt2/2, -b):
accB=Acc(); int_simple_kernel(accB, fneg(hr2), smb, P12)
dB = mp.quad(lambda y: -(r2/2)/(y+b)*hg1(y)*hg2(y), [0,1])
print("kernel -(r2/2)/(y+b): engine", mp.nstr(acc_value(accB),18), " direct", mp.nstr(dB,18))
