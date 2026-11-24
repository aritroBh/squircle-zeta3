"""RST1 step_05: precise frequency/phase of the rest oscillation.
Compare arg(R_k) with theta2=sqrt(2 pi N) and with the saddle phase ImPhi(s_*).
Local frequency check: arg(R_{k+1})-arg(R_k) vs sqrt(2pi)(sqrt(N+1)-sqrt(N)).
"""
import mpmath as mp
mp.mp.dps = 80
pi=mp.pi; j=1j; omega=mp.expjpi(mp.mpf(1)/4)
def H(s): return mp.beta(s/2, mp.mpf(3)/2)/2
sq2pi=mp.sqrt(2*pi)

def Phi(s,k):
    return j*pi*s/2 + (k-1)*mp.log(s) + mp.log(s+1) - (k+1)*mp.log(s+2)
def dPhi(s,k):
    return j*pi/2 + (k-1)/s + 1/(s+1) - (k+1)/(s+2)

def Rmodel(k):
    rs=2*mp.sqrt(k/pi)
    def f(r):
        s=r*omega
        return H(s)*mp.expj(pi*s/2)*s**(k-1)*(s+1)/(s+2)**(k+1)*omega
    pts=[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf]
    return -4/pi*mp.quad(f,pts).imag

# unwrapped phase of R_k
prev=None; un=None
print(" k     arg(R_k) unwrapped - theta2      ImPhi(s*) - theta2 (mod 2pi, centered)")
for k in range(380, 421, 4):
    R=Rmodel(k); N=k+1; th=sq2pi*mp.sqrt(N)
    a=mp.arg(R)
    if un is None: un=a
    else:
        da=a-prev
        if da>pi: da-=2*pi
        if da<-pi: da+=2*pi
        un+=da
    prev=a
    r0=2*mp.sqrt((k-1)/pi)
    s_=mp.findroot(lambda s: dPhi(s,k), r0*omega)
    imP = mp.fmod(mp.im(Phi(s_,k))-th+pi, 2*pi)-pi
    print(f"{k:4d}   {mp.nstr(un-th,10):>14s}      {mp.nstr(imP,10)}")
# local frequency
print()
print("local frequency: (un(k+4)-un(k))/4 vs (theta2(k+4)-theta2(k))/4")
