import sys, json
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
def Iplus(k,a,jj):
    N=k+1; E=mp.e**(1j*a); gp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))
    g1=-mp.mpf(3)/4-1j*E/(4*mp.sin(a)); gj=[mp.mpc(1),g1][jj]
    xs=mp.sqrt(N*a)
    J=mp.quad(lambda x:mp.e**(-x)*x**(N-1)*(x+1j*a)**(-(N-jj-mp.mpf(3)/2)),[0,xs/4,xs,4*xs,mp.inf])
    return gp*gj*(1j/mp.pi)*mp.gamma(jj+mp.mpf(3)/2)*mp.gamma(N-jj-mp.mpf(3)/2)/mp.gamma(N)*J
a=2*mp.pi/3
M=[mp.mpf(s) for s in json.load(open('mk_2pi3_k80.json'))]
print("k : |M_k - 2Re(I0+I1)| / (N^-2 e^-theta)   [Watson remainder constant]")
for k in [40,50,60,70,80]:
    N=k+1; th=mp.sqrt(2*a*N)
    main=2*mp.re(Iplus(k,a,0)+Iplus(k,a,1))
    rest=abs(M[k]-main)
    print(f"k={k}: rest={mp.nstr(rest,4)}  /(N^-2 e^-th)={mp.nstr(rest/(N**(-2)*mp.e**(-th)),5)}  /(N^-3/2 e^-th)=sc {mp.nstr(rest/(N**(-1.5)*mp.e**(-th)),5)}")
