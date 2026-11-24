import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 55
from lead2r_lib import Efun
h = lambda w: ((1+w)/mp.sqrt(1+w**2)-1)/w
q = lambda w: (1+w)/(1+w**2)**mp.mpf('1.5')
g1 = lambda w: (1-w)/(1+w)*h(w)
cn = lambda n: mp.binomial(2*n,n)/4**n
def F_series(s):  # F(s) = int_0^1 w^{-s} h dw via proved partial fractions
    return mp.nsum(lambda n: (-1)**int(n)*cn(n)/(2*n-s), [1, mp.inf]) + \
           mp.nsum(lambda n: (-1)**int(n)*cn(n)/(2*n+1-s), [0, mp.inf])
def int01_sing(s, f):
    """int_0^1 x^s f(x) dx, f(0)=1, Re s > -1: geometric slices + explicit tail eps^{s+1}/(s+1)."""
    tot = mp.mpf(0); x = mp.mpf(1)
    while x > mp.mpf(10)**(-40):
        tot += mp.quad(lambda t: t**s*f(t), [x/2, x]); x /= 2
    return tot + x**(s+1)/(s+1)
pts = [mp.mpf('0.3'), mp.mpf('-0.4'), mp.mpf('0.73'), mp.mpf('-0.9'), mp.mpf('0.99'), mp.mpf('-0.99'),
       mp.mpc('0.2','0.4'), mp.mpc('-0.35','-0.28'), mp.mpc('0.01','0.9'), mp.mpc('-0.97','0.15')]
worst = {k: mp.mpf(0) for k in ('WW','FF','FFser','GG','IBP')}
for s in pts:
    WW = int01_sing(-s,q) + int01_sing(s,q);   rW = Efun(s)+Efun(-s)
    F1 = int01_sing(-s,h); F2 = int01_sing(s,h)
    FF = F1+F2; rF = (Efun(s)-Efun(-s))/s
    GG = int01_sing(-s,g1) - int01_sing(s,g1); rG = 2*mp.pi/mp.sin(mp.pi*s)-(Efun(s)+Efun(-s))/s
    IBP = (Efun(s)-Efun(-s)) - s*FF
    Fs = F_series(s)
    worst['WW']=max(worst['WW'],abs(WW-rW)); worst['FF']=max(worst['FF'],abs(FF-rF))
    worst['GG']=max(worst['GG'],abs(GG-rG)); worst['IBP']=max(worst['IBP'],abs(IBP))
    worst['FFser']=max(worst['FFser'],abs(F1-Fs))
for k, val in worst.items(): print(f"worst {k}: {mp.nstr(val,3)}")
with open('progress.log','a') as f:
    f.write("step_13: Master Lemma all parts at 10 strip pts (incl +-0.99, complex) dps55, slices+tail quadrature AND partial-fraction series: WW/FF/GG/IBP worst "+mp.nstr(max(worst['WW'],worst['FF'],worst['GG'],worst['IBP']),3)+"; F-series==quad worst "+mp.nstr(worst['FFser'],3)+"\n")
