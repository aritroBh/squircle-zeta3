import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 50
from lead2r_lib import Efun
h = lambda w: ((1+w)/mp.sqrt(1+w**2)-1)/w
g1 = lambda w: (1-w)/(1+w)*h(w)
def J(t, f, R=None):
    # split at R = c/|t+1| to control slow exponential decay
    a = abs(t+1)
    if R is None:
        R = mp.mpf(30)/a if a > mp.mpf('1e-30') else mp.mpf(100)
    return mp.quad(lambda r: mp.exp(-(t+1)*r)*f(mp.exp(-r)), [0, R]) + \
           mp.quad(lambda r: mp.exp(-(t+1)*r)*f(mp.exp(-r)), [R, mp.inf])
pts = [mp.mpf('0.3'), mp.mpf('-0.4'), mp.mpf('0.73'), mp.mpf('-0.9'), mp.mpf('0.99'), mp.mpf('-0.99'),
       mp.mpc('0.2','0.4'), mp.mpc('-0.35','-0.28'), mp.mpc('0.01','0.9'), mp.mpc('-0.97','0.15')]
worst = {k: mp.mpf(0) for k in ('FF','GG','IBP')}
for s in pts:
    FF = J(-s,h) + J(s,h);  rF = (Efun(s)-Efun(-s))/s
    GG = J(-s,g1) - J(s,g1); rG = 2*mp.pi/mp.sin(mp.pi*s)-(Efun(s)+Efun(-s))/s
    IBP = (Efun(s)-Efun(-s)) - s*FF
    for k,v in (('FF',abs(FF-rF)),('GG',abs(GG-rG)),('IBP',abs(IBP))):
        worst[k] = max(worst[k], v)
        if v > mp.mpf('1e-30'): print(f"  s={mp.nstr(s,8)} {k}: {mp.nstr(v,3)}")
for k, val in worst.items(): print(f"worst {k}: {mp.nstr(val,3)}")
with open('progress.log','a') as f:
    f.write("step_11: FF/GG/IBP with decay-adaptive split: worst "+mp.nstr(max(worst.values()),3)+"\n")
