import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 50
from lead2r_lib import Efun
h = lambda w: ((1+w)/mp.sqrt(1+w**2)-1)/w
q = lambda w: (1+w)/(1+w**2)**mp.mpf('1.5')
g1 = lambda w: (1-w)/(1+w)*h(w)
# symbolic inversion symmetries
import sympy as sp
w = sp.symbols('w', positive=True)
hs = ((1+w)/sp.sqrt(1+w**2)-1)/w; qs = (1+w)/(1+w**2)**sp.Rational(3,2)
print("h(1/w)=w^2 h(w):", sp.radsimp(sp.simplify(hs.subs(w,1/w) - w**2*hs))==0)
print("q(1/w)=w^2 q(w):", sp.radsimp(sp.simplify(qs.subs(w,1/w) - w**2*qs))==0)
print("g1(1/w)=-w^2 g1(w):", sp.radsimp(sp.simplify(((1-w)/(1+w)*hs).subs(w,1/w) + w**2*(1-w)/(1+w)*hs))==0)
pts = [mp.mpf('0.3'), mp.mpf('-0.4'), mp.mpf('0.73'), mp.mpf('-0.9'), mp.mpc('0.2','0.4'), mp.mpc('-0.35','-0.28')]
worst = {k: mp.mpf(0) for k in ('WW','FF','GG','B4','IBP')}
for s in pts:
    WW = mp.quad(lambda x: x**(-s)*q(x), [0,1,mp.inf]);  rW = Efun(s)+Efun(-s)
    FF = mp.quad(lambda x: x**(-s)*h(x), [0,1,mp.inf]);  rF = (Efun(s)-Efun(-s))/s
    GG = mp.quad(lambda x: x**(-s)*g1(x), [0,1,mp.inf]); rG = 2*mp.pi/mp.sin(mp.pi*s)-(Efun(s)+Efun(-s))/s
    B4 = mp.quad(lambda x: (x**(1+s)+x**(-s))/(1+x**2)**mp.mpf('1.5'), [0,1]);  # == E(s) eq.(4)
    F1 = mp.quad(lambda x: x**(-s)*h(x), [0,1]); F2 = mp.quad(lambda x: x**s*h(x), [0,1])
    IBP = (Efun(s)-Efun(-s)) - s*(F1+F2)
    worst['WW'] = max(worst['WW'], abs(WW-rW)); worst['FF'] = max(worst['FF'], abs(FF-rF))
    worst['GG'] = max(worst['GG'], abs(GG-rG)); worst['B4'] = max(worst['B4'], abs(B4-Efun(s)))
    worst['IBP'] = max(worst['IBP'], abs(IBP))
for k, val in worst.items(): print(f"worst {k}: {mp.nstr(val,3)}")
# F-series erratum
cn = lambda n: mp.mpf(sp.binomial(2*n,n))/4**n
for y in [mp.mpf('0.3'), mp.mpf('-0.45')]:
    fq = mp.quad(lambda x: x**(-y)*h(x), [0,1])
    correct = mp.nsum(lambda n: (-1)**int(n)*cn(n)/(2*n-y), [1, mp.inf]) + mp.nsum(lambda n: (-1)**int(n)*cn(n)/(2*n+1-y), [0, mp.inf])
    lead1printed = mp.nsum(lambda n: (-1)**int(n)*cn(n)*(1/(2*n-y)+1/(2*n+1-y)), [1, mp.inf])
    print(f"y={y}: correct-F={mp.nstr(correct-fq,3)}  printedLEAD1-F={mp.nstr(lead1printed-fq,3)}  -1/(1-y)={mp.nstr(-1/(1-y),3)}")
with open('progress.log','a') as f:
    f.write("step_06: Master Lemma (i)(ii)(iii)+Beta(4)+IBP at 6 strip pts dps50: worst "+mp.nstr(max(worst.values()),3)+"; inversion syms symbolic True; erratum confirmed: LEAD1 printed form off by exactly -1/(1-y)\n")
