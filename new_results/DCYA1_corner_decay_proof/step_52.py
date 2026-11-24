import sys, json
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
def Iplus(k, a, order=0):
    N=k+1; E=mp.e**(1j*a)
    gp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))
    g1=-mp.mpf(3)/4-1j*E/(4*mp.sin(a))
    xs=mp.sqrt(N*a)
    def Jj(jj):  # int e^-x x^{N-1} (x+i a)^{-(N-jj-3/2)} dx
        f=lambda x: mp.e**(-x)*x**(N-1)*(x+1j*a)**(-(N-jj-mp.mpf(3)/2))
        return mp.quad(f,[0,xs/4,xs,4*xs,mp.inf])
    I0=gp*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)*Jj(0)
    if order==0: return I0
    I1=gp*g1*(1j/mp.pi)*mp.gamma(mp.mpf(5)/2)*mp.gamma(N-mp.mpf(5)/2)/mp.gamma(N)*Jj(1)
    return I0+I1
for tag,a in [('2pi3',2*mp.pi/3),('pi3',mp.pi/3)]:
    M=[mp.mpf(s) for s in json.load(open(f'mk_{tag}_k80.json'))]
    print(f"=== alpha={tag} ===")
    for k in [60,80]:
        m0=2*mp.re(Iplus(k,a,0)); m1=2*mp.re(Iplus(k,a,1)); N=k+1; th=mp.sqrt(2*a*N); sc=N**(-1.5)*mp.e**(-th)
        print(f" k={k}: |Mk-main0|/sc={mp.nstr(abs(M[k]-m0)/sc,4)}  |Mk-main0-g1|/sc={mp.nstr(abs(M[k]-m1)/sc,4)}")
