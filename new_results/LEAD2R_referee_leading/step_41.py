import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 50
from lead2r_lib import Efun
h = lambda w: ((1+w)/mp.sqrt(1+w**2)-1)/w
g1 = lambda w: (1-w)/(1+w)*h(w)
def Jr(t, f):
    # int_0^inf e^{-(t+1)r} f(e^{-r}) dr = 1/(t+1) + int e^{-(t+1)r} (f(e^{-r})-1) dr ; remainder decays e^{-(t+2)r}
    return 1/(t+1) + mp.quad(lambda r: mp.exp(-(t+1)*r)*(f(mp.exp(-r))-1), [0, 20, mp.inf])
pts = [mp.mpf('0.3'), mp.mpf('-0.4'), mp.mpf('0.73'), mp.mpf('-0.9'), mp.mpf('0.99'), mp.mpf('-0.99'),
       mp.mpc('0.2','0.4'), mp.mpc('-0.35','-0.28'), mp.mpc('0.01','0.9'), mp.mpc('-0.97','0.15')]
worst = {k: mp.mpf(0) for k in ('FF','GG','IBP')}
for s in pts:
    FF = Jr(-s,h) + Jr(s,h);  rF = (Efun(s)-Efun(-s))/s
    GG = Jr(-s,g1) - Jr(s,g1); rG = 2*mp.pi/mp.sin(mp.pi*s)-(Efun(s)+Efun(-s))/s
    IBP = (Efun(s)-Efun(-s)) - s*FF
    worst['FF']=max(worst['FF'],abs(FF-rF)); worst['GG']=max(worst['GG'],abs(GG-rG)); worst['IBP']=max(worst['IBP'],abs(IBP))
for k, val in worst.items(): print(f"worst {k}: {mp.nstr(val,3)}")
with open('progress.log','a') as f:
    f.write("step_12: FF/GG/IBP with asymptotic subtraction: all 10 strip pts dps50 worst "+mp.nstr(max(worst.values()),3)+"\n")
