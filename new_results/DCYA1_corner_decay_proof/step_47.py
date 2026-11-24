import time, sys, json
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
def H(s): return mp.beta(s/2, mp.mpf(3)/2)/2
def Mmain_pi2(k):
    w=mp.e**(1j*mp.pi/4)
    def f(r):
        s=r*w
        Q=(s/(1+s))**k/(1+s)
        return Q*H(s)*mp.e**(1j*mp.pi*s/2)*w
    rs=mp.sqrt(2*(k+1)/mp.pi)
    I=mp.quad(f,[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf])
    return (2/mp.pi)*mp.im(I)   # + sign
M=[mp.mpf(s) for s in json.load(open('mk_pi2_k100_dps140.json'))]
for k in [40,60,80,100]:
    mm=Mmain_pi2(k); N=k+1; th=mp.sqrt(mp.pi*N)
    lead=-(mp.sqrt(mp.pi)/N)*mp.e**(-th)*mp.cos(th)
    sc=N**(-1.5)*mp.e**(-th)
    print(f"k={k}: Mmain={mp.nstr(mm,13)} Mk={mp.nstr(M[k],13)} |Mk-Mmain|/sc={mp.nstr(abs(M[k]-mm)/sc,4)} |Mmain-lead|/sc={mp.nstr(abs(mm-lead)/sc,4)}")
