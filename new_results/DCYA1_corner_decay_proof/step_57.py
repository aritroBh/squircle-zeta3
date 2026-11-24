import sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=40
a=2*mp.pi/3
def J0(N,a):
    f=lambda x: mp.e**(-x)*x**(N-1)*(x+1j*a)**(-(N-mp.mpf(3)/2))
    xs=mp.sqrt(N*a)
    return mp.quad(f,[0,xs/4,xs,4*xs,mp.inf])
prev=None
for N in [200,400,800,1600,3200]:
    th=mp.sqrt(2*a*N)
    v=J0(N,a)
    env=abs(v)*mp.e**th  # ~ N^p
    print(f"N={N}: |J0|e^th={mp.nstr(env,6)}  N^(1/2)={mp.nstr(mp.sqrt(N),4)}  ratio env/sqrt(N)={mp.nstr(env/mp.sqrt(N),5)}")
