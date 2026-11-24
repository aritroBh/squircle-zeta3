import sys, json
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
def H(s): return mp.beta(s/2, mp.mpf(3)/2)/2
def Mmain(k, a, Ha):
    w=mp.e**(1j*mp.pi/4)
    def f(r):
        s=r*w
        Q=(s/(1+s))**k/(1+s)
        return Q*Ha(s)*mp.e**(1j*a*s)*w
    rs=mp.sqrt(2*(k+1)/a)
    I=mp.quad(f,[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf])
    return (2/mp.pi)*mp.im(I)
for tag,a in [('2pi3',2*mp.pi/3),('pi3',mp.pi/3)]:
    M=[mp.mpf(s) for s in json.load(open(f'mk_{tag}_k80.json'))]
    A=mp.sqrt(2*a*mp.sin(a))
    print(f"=== alpha={tag} ===")
    for k in [40,60,80]:
        mm=Mmain(k,a,H); N=k+1; th=mp.sqrt(2*a*N); sc=N**(-1.5)*mp.e**(-th)
        lead=-(A/N)*mp.e**(-th)*mp.cos(th)
        print(f" k={k}: Mmain(H)={mp.nstr(mm,10)} Mk={mp.nstr(M[k],10)} Mk/Mmain={mp.nstr(M[k]/mm,6)} |Mmain-lead|/sc={mp.nstr(abs(mm-lead)/sc,3)}")
