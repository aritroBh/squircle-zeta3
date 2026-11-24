import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 50
from lead2r_lib import Efun
h = lambda w: ((1+w)/mp.sqrt(1+w**2)-1)/w
q = lambda w: (1+w)/(1+w**2)**mp.mpf('1.5')
g1 = lambda w: (1-w)/(1+w)*h(w)
pts = [mp.mpf('0.3'), mp.mpf('-0.4'), mp.mpf('0.73'), mp.mpf('-0.9'), mp.mpf('0.99'), mp.mpc('0.2','0.4'), mp.mpc('-0.35','-0.28'), mp.mpc('0.01','0.9')]
for s in pts:
    WW = mp.quad(lambda x: (x**s+x**(-s))*q(x), [0,1]);  rW = Efun(s)+Efun(-s)
    FF = mp.quad(lambda x: (x**s+x**(-s))*h(x), [0,1]);  rF = (Efun(s)-Efun(-s))/s
    GG = mp.quad(lambda x: (x**(-s)-x**s)*g1(x), [0,1]); rG = 2*mp.pi/mp.sin(mp.pi*s)-(Efun(s)+Efun(-s))/s
    print(f"s={mp.nstr(s,10)}: WW {mp.nstr(abs(WW-rW),2)} FF {mp.nstr(abs(FF-rF),2)} GG {mp.nstr(abs(GG-rG),2)}")
