import os; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log as flog
mp.mp.dps = 100
r2=mp.sqrt(2); al=r2-1; be=r2/2-1; ln=mp.log
def Li21(x,y):
    x=mp.mpf(x); y=mp.mpf(y); N=int(mp.mp.dps/max(mp.log10(1/abs(x)),mp.mpf('0.05')))+60
    s=mp.mpf(0); Hy=mp.mpf(0); py=mp.mpf(1)
    for n in range(1,N+1):
        s += x**n/n**2*Hy; py*=y; Hy+=py/n
    return s
U = mp.quad(lambda t: ln(t)*ln(1+be*t)/(1+al*t), [mp.mpf(0),mp.mpf('0.5'),1])
pred = -(1/al)*Li21(-al, be/al)
print("U        =", mp.nstr(U,40))
print("-(1/a)Li21(-a,b/a) =", mp.nstr(pred,40))
print("diff =", mp.nstr(abs(U-pred),5))
# also -al = x1, be/al = ?
print("-al =", mp.nstr(-al,10), "= 1-sqrt2?", mp.nstr(1-r2,10))
print("be/al =", mp.nstr(be/al,10), "= -sqrt2/2?", mp.nstr(-r2/2,10))
flog("step: verified U = -(1/a)Li21(x1,x2/x1), diff=%s"%mp.nstr(abs(U-pred),3))
