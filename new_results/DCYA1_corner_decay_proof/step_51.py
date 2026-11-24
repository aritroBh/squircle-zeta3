import sys, json
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
# K2-derived 1D main term: M_main = 2 Re[ gamma_+ * (i/pi) * Gamma(3/2) * Gamma(N-3/2)/Gamma(N) * J ]
# J = int_0^inf e^{-x} x^{N-1} (x+i a)^{-(N-3/2)} dx
def Iplus(k, a):
    N=k+1
    gp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))
    def f(x):
        return mp.e**(-x)*x**(N-1)*(x+1j*a)**(-(N-mp.mpf(3)/2))
    xs=mp.sqrt(N*a)
    J=mp.quad(f,[0,xs/4,xs,4*xs,mp.inf])
    return gp*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)*J
for tag,a in [('2pi3',2*mp.pi/3),('pi3',mp.pi/3)]:
    M=[mp.mpf(s) for s in json.load(open(f'mk_{tag}_k80.json'))]
    A=mp.sqrt(2*a*mp.sin(a))
    print(f"=== alpha={tag} ===")
    for k in [40,60,80]:
        mm=2*mp.re(Iplus(k,a)); N=k+1; th=mp.sqrt(2*a*N); sc=N**(-1.5)*mp.e**(-th)
        print(f" k={k}: Mmain_K2={mp.nstr(mm,11)} Mk={mp.nstr(M[k],11)} |Mk-Mmain|/sc={mp.nstr(abs(M[k]-mm)/sc,4)}")
