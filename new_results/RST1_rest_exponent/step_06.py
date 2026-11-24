"""RST1 step_06: large-k saddle constants for the rest asymptotics.
amplitude: |GM|*N e^{theta2};  phase: ImPhi(s*)-theta2;  total constant via gk.
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
def ddPhi(s,k):
    return -(k-1)/s**2 - 1/(s+1)**2 + (k+1)/(s+2)**2
print(" k        |GM|*N e^th2        ImPhi-th2          arg(gkf) - th2 (mod2pi, centered)")
for k in [400,1600,6400,25600,102400,409600]:
    N=k+1; th=sq2pi*mp.sqrt(N)
    r0=2*mp.sqrt((k-1)/pi)
    s_=mp.findroot(lambda s: dPhi(s,k), r0*omega)
    g = omega*mp.exp(Phi(s_,k))*H(s_)*mp.sqrt(2*pi/(-ddPhi(s_,k)))
    GM = -(4/pi)*g.imag
    ph = mp.fmod(mp.im(Phi(s_,k))-th+pi,2*pi)-pi
    # gkf = -(4/pi)*g (complex before Im); its argument determines cosine phase
    gkf = -(4/pi)*g
    ag = mp.fmod(mp.arg(gkf)-th+pi,2*pi)-pi
    print(f"{k:7d}  {mp.nstr(abs(GM)*N*mp.e**th,10):>14s}   {mp.nstr(ph,12):>14s}   {mp.nstr(ag,12)}")
print()
print("candidates: sqrt(pi)=",mp.nstr(mp.sqrt(pi),8)," 4/sqrt(pi)=",mp.nstr(4/mp.sqrt(pi),8))
print("phase candidates: -3pi/4=",mp.nstr(-3*pi/4,8)," -pi+1/2=",mp.nstr(-pi+mp.mpf('0.5'),8))
