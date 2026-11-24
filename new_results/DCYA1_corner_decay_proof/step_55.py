import sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
def Iplus0(k, a):
    N=k+1
    gp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))
    f=lambda x: mp.e**(-x)*x**(N-1)*(x+1j*a)**(-(N-mp.mpf(3)/2))
    xs=mp.sqrt(N*a)
    J=mp.quad(f,[0,xs/4,xs,4*xs,mp.inf])
    return gp*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)*J
a=2*mp.pi/3; A=mp.sqrt(2*a*mp.sin(a))
print("N-scaling of 2 Re I_+^{(0)} and of lead:")
for k in [100,200,400,800]:
    N=k+1; th=mp.sqrt(2*a*N)
    v=2*mp.re(Iplus0(k,a))
    lead=-(A/N)*mp.e**(-th)*mp.cos(th)
    # fit power: |v| e^{th} ~ N^p (use envelope, divide by |lead|*N/A = e^-th)
    env=abs(v)*mp.e**th
    print(f" k={k}: 2ReI0={mp.nstr(v,8)}  lead={mp.nstr(lead,8)}  ratio v/lead={mp.nstr(v/lead,5)}  |v|e^th={mp.nstr(env,5)}")
