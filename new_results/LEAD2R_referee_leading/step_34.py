import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 50
from lead2r_lib import Efun
h = lambda w: ((1+w)/mp.sqrt(1+w**2)-1)/w
q = lambda w: (1+w)/(1+w**2)**mp.mpf('1.5')
g1 = lambda w: (1-w)/(1+w)*h(w)
# reduced proper forms (derived from the proved symmetries): WW=int_0^1 (x^s+x^-s) q, FF=int_0^1 (x^s+x^-s) h, GG=int_0^1 (x^-s - x^s) g1
pts = [mp.mpf('0.3'), mp.mpf('-0.4'), mp.mpf('0.73'), mp.mpf('-0.9'), mp.mpf('0.99'), mp.mpc('0.2','0.4'), mp.mpc('-0.35','-0.28'), mp.mpc('0.01','0.9')]
worst = {k: mp.mpf(0) for k in ('WW','FF','GG','GGdirect')}
for s in pts:
    WW = mp.quad(lambda x: (x**s+x**(-s))*q(x), [0,1]);  rW = Efun(s)+Efun(-s)
    FF = mp.quad(lambda x: (x**s+x**(-s))*h(x), [0,1]);  rF = (Efun(s)-Efun(-s))/s
    GG = mp.quad(lambda x: (x**(-s)-x**s)*g1(x), [0,1]); rG = 2*mp.pi/mp.sin(mp.pi*s)-(Efun(s)+Efun(-s))/s
    # GG also directly on [0,inf] with split for cross-check
    GGd = mp.quad(lambda x: x**(-s)*g1(x), [0,1]) + mp.quad(lambda x: x**(-s)*g1(x), [1, mp.inf])
    worst['WW'] = max(worst['WW'], abs(WW-rW)); worst['FF'] = max(worst['FF'], abs(FF-rF))
    worst['GG'] = max(worst['GG'], abs(GG-rG)); worst['GGdirect'] = max(worst['GGdirect'], abs(GGd-rG))
for k, val in worst.items(): print(f"worst {k}: {mp.nstr(val,3)}")
with open('progress.log','a') as f:
    f.write("step_07: Master Lemma via proper [0,1] forms at 8 strip pts dps50: WW/FF/GG worst "+mp.nstr(max(worst['WW'],worst['FF'],worst['GG']),3)+" (direct inf-form GG worst "+mp.nstr(worst['GGdirect'],3)+")\n")
