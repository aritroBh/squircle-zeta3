import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 50
from lead2r_lib import Efun
h = lambda w: ((1+w)/mp.sqrt(1+w**2)-1)/w
q = lambda w: (1+w)/(1+w**2)**mp.mpf('1.5')
g1 = lambda w: (1-w)/(1+w)*h(w)
qb = lambda w: (1+w**2)**mp.mpf('-1.5')
def I(s, f):   # int_0^1 x^{-s} f(x) dx, f(0)=1, regularized
    return 1/(1-s) + mp.quad(lambda x: x**(-s)*(f(x)-1), [0,1])
pts = [mp.mpf('0.3'), mp.mpf('-0.4'), mp.mpf('0.73'), mp.mpf('-0.9'), mp.mpf('0.99'), mp.mpf('-0.99'),
       mp.mpc('0.2','0.4'), mp.mpc('-0.35','-0.28'), mp.mpc('0.01','0.9'), mp.mpc('-0.97','0.15')]
worst = {k: mp.mpf(0) for k in ('FF','GG','B4','IBP')}
for s in pts:
    F1 = I(s,h); F2 = I(-s,h)   # note: int_0^1 x^{s} h = I(-s, h) regularized the other way
    FF = F1+F2; rF = (Efun(s)-Efun(-s))/s
    G1 = I(s,g1); G2 = I(-s,g1)
    GG = G1-G2; rG = 2*mp.pi/mp.sin(mp.pi*s)-(Efun(s)+Efun(-s))/s
    B4 = I(-1-s, qb) if False else (I(s, lambda x: x*qb(x)) + I(s, qb))  # split: x^{1+s}qb + x^{-s}qb
    IBP = (Efun(s)-Efun(-s)) - s*(F1+F2)
    for k,v in (('FF',abs(FF-rF)),('GG',abs(GG-rG)),('B4',abs(B4-Efun(s))),('IBP',abs(IBP))):
        if v > worst[k]:
            worst[k]=v
            print(f"  new worst {k} at s={mp.nstr(s,8)}: {mp.nstr(v,3)}")
for k, val in worst.items(): print(f"worst {k}: {mp.nstr(val,3)}")
with open('progress.log','a') as f:
    f.write("step_09: Master Lemma FF/GG/B4/IBP fully regularized at 10 pts dps50: worst "+mp.nstr(max(worst.values()),3)+"\n")
