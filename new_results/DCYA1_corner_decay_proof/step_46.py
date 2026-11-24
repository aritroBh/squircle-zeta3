import time, sys, json
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
# D13-style main term at pi/2: M_main = -(2/pi) Im int_0^inf Q_k(s) H(s) e^{i pi s/2} ds, ray s=e^{i pi/4} x
def H(s): return mp.beta(s/2, mp.mpf(3)/2)/2
def Mmain_pi2(k):
    w=mp.e**(1j*mp.pi/4)
    def f(r):
        s=r*w
        Q=(s/(1+s))**k/(1+s)   # = s^k/(1+s)^{k+1}
        return Q*H(s)*mp.e**(1j*mp.pi*s/2)*w
    rs=mp.sqrt(2*(k+1)/mp.pi)
    I=mp.quad(f,[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf])
    return -(2/mp.pi)*mp.im(I)
# accurate M_k pi/2 (dps=140 file)
M=[mp.mpf(s) for s in json.load(open('mk_pi2_k100_dps140.json'))]
for k in [40,60,80,100]:
    mm=Mmain_pi2(k)
    N=k+1; th=mp.sqrt(mp.pi*N); lead=-(mp.sqrt(mp.pi)/N)*mp.e**(-th)*mp.cos(th)
    print(f"k={k}: Mmain={mp.nstr(mm,14)}")
    print(f"      Mk   ={mp.nstr(M[k],14)}  Mk-Mmain={mp.nstr(abs(M[k]-mm),3)}")
    print(f"      lead ={mp.nstr(lead,14)}  |Mmain-lead|/(N^-3/2 e^-th)={mp.nstr(abs(mm-lead)/(N**(-1.5)*mp.e**(-th)),4)}")
