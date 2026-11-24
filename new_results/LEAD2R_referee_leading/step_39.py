import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 50
from lead2r_lib import Efun
h = lambda w: ((1+w)/mp.sqrt(1+w**2)-1)/w
q = lambda w: (1+w)/(1+w**2)**mp.mpf('1.5')
g1 = lambda w: (1-w)/(1+w)*h(w)
qb = lambda w: (1+w**2)**mp.mpf('-1.5')
# int_0^1 x^{t} f(x) dx = int_0^inf e^{-(t+1)r} f(e^{-r}) dr  (Re t > -1)
def J(t, f):
    return mp.quad(lambda r: mp.exp(-(t+1)*r)*f(mp.exp(-r)), [0, mp.inf])
pts = [mp.mpf('0.3'), mp.mpf('-0.4'), mp.mpf('0.73'), mp.mpf('-0.9'), mp.mpf('0.99'), mp.mpf('-0.99'),
       mp.mpc('0.2','0.4'), mp.mpc('-0.35','-0.28'), mp.mpc('0.01','0.9'), mp.mpc('-0.97','0.15')]
worst = {k: mp.mpf(0) for k in ('WW','FF','GG','B4','IBP')}
for s in pts:
    WW = J(-s,q) + J(s,q);              rW = Efun(s)+Efun(-s)
    FF = J(-s,h) + J(s,h);              rF = (Efun(s)-Efun(-s))/s
    GG = J(-s,g1) - J(s,g1);            rG = 2*mp.pi/mp.sin(mp.pi*s)-(Efun(s)+Efun(-s))/s
    B4 = J(1+s,qb) + J(-s,qb)           # eq (4): E(s) = int_0^1 (w^{1+s}+w^{-s})/(1+w^2)^{3/2} dw
    IBP = (Efun(s)-Efun(-s)) - s*FF
    worst['WW']=max(worst['WW'],abs(WW-rW)); worst['FF']=max(worst['FF'],abs(FF-rF))
    worst['GG']=max(worst['GG'],abs(GG-rG)); worst['B4']=max(worst['B4'],abs(B4-Efun(s)))
    worst['IBP']=max(worst['IBP'],abs(IBP))
for k, val in worst.items(): print(f"worst {k}: {mp.nstr(val,3)}")
with open('progress.log','a') as f:
    f.write("step_10: Master Lemma via x=e^-r smooth quadrature, 10 strip pts incl +-0.99, dps50: WW/FF/GG/B4/IBP worst "+mp.nstr(max(worst.values()),3)+"\n")
