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
tag='2pi3'; a=2*mp.pi/3
M=[mp.mpf(s) for s in json.load(open(f'mk_{tag}_k80.json'))]
A=mp.sqrt(2*a*mp.sin(a))
cands={
 'H': lambda s: H(s),
 'sqrt(2a sina/pi) H': lambda s: mp.sqrt(2*a*mp.sin(a)/mp.pi)*H(s),
 'sqrt(sina) H': lambda s: mp.sqrt(mp.sin(a))*H(s),
 'H(a s/(pi/2))': lambda s: H(a*s/(mp.pi/2)),
}
for nm,Ha in cands.items():
    r=[M[k]/Mmain(k,a,Ha) for k in [60,80]]
    print(f"{nm:26s}: Mk/Mmain k=60,80 = {[mp.nstr(x,5) for x in r]}")
