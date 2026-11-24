import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 50
from lead2r_lib import Efun
h = lambda w: ((1+w)/mp.sqrt(1+w**2)-1)/w
q = lambda w: (1+w)/(1+w**2)**mp.mpf('1.5')
g1 = lambda w: (1-w)/(1+w)*h(w)
def int_xs_f(s, f):
    # int_0^1 x^{-s} f(x) dx with f(0)=1: subtract the singular part
    return 1/(1-s) + mp.quad(lambda x: x**(-s)*(f(x)-1), [0,1])
def int_xs_f_plain(s, f):  # int_0^1 x^{s} f dx (Re s > -1, no singularity issue unless near -1)
    return 1/(1+s) + mp.quad(lambda x: x**(s)*(f(x)-1), [0,1])
pts = [mp.mpf('0.3'), mp.mpf('-0.4'), mp.mpf('0.73'), mp.mpf('-0.9'), mp.mpf('0.99'), mp.mpf('-0.99'),
       mp.mpc('0.2','0.4'), mp.mpc('-0.35','-0.28'), mp.mpc('0.01','0.9'), mp.mpc('-0.97','0.15')]
worst = {k: mp.mpf(0) for k in ('WW','FF','GG','B4','IBP','Mgpieces')}
for s in pts:
    WW = int_xs_f(s,q) + int_xs_f_plain(s,q);  rW = Efun(s)+Efun(-s)
    FF = int_xs_f(s,h) + int_xs_f_plain(s,h);  rF = (Efun(s)-Efun(-s))/s
    GG = int_xs_f(s,g1) - int_xs_f_plain(s,g1); rG = 2*mp.pi/mp.sin(mp.pi*s)-(Efun(s)+Efun(-s))/s
    B4 = mp.quad(lambda x: (x**(1+s)+x**(-s))/(1+x**2)**mp.mpf('1.5'), [0,1])
    F1 = int_xs_f(s,h); F2 = int_xs_f_plain(s,h)
    IBP = (Efun(s)-Efun(-s)) - s*(F1+F2)
    # Mellin pieces for GG: M_g(t) = Gamma(t/2)Gamma((1-t)/2)/(2sqrt pi) - pi/sin(pi t), and GG(s)=M_g(-s)-M_g(1-s)
    t = -s
    Mg1 = mp.quad(lambda x: x**(t-1)/mp.sqrt(1+x**2), [0,1,mp.inf]) if abs(mp.re(t))<0.9 else None
    worst['WW'] = max(worst['WW'], abs(WW-rW)); worst['FF'] = max(worst['FF'], abs(FF-rF))
    worst['GG'] = max(worst['GG'], abs(GG-rG)); worst['B4'] = max(worst['B4'], abs(B4-Efun(s)))
    worst['IBP'] = max(worst['IBP'], abs(IBP))
for k, val in worst.items(): print(f"worst {k}: {mp.nstr(val,3)}")
with open('progress.log','a') as f:
    f.write("step_08: Master Lemma (i)(ii)(iii)+Beta(4)+IBP at 10 strip pts incl +-0.99, complex, dps50 with regularized quadrature: worst "+mp.nstr(max(worst.values()),3)+"\n")
