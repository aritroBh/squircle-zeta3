"""RST1 step_04: saddle asymptotics of the q=2 rest model R_k.
R_k = -(4/pi) Im int_ray H e^{i pi s/2} s^{k-1}(s+1)/(s+2)^{k+1} ds.
Saddle s_* solves Phi'=0; Gaussian model GM = -(4/pi) Im[omega e^{Phi(s_*)} H(s_*) sqrt(2pi/(-Phi''(s_*)))].
Fit leading form: R_k ~ C/N e^{-sqrt(2 pi N)} cos(sqrt(2 pi N) + phi).
"""
import mpmath as mp
mp.mp.dps = 90
pi=mp.pi; j=1j; omega=mp.expjpi(mp.mpf(1)/4)
def H(s): return mp.beta(s/2, mp.mpf(3)/2)/2
sq2pi=mp.sqrt(2*pi)

def Phi(s,k):
    return j*pi*s/2 + (k-1)*mp.log(s) + mp.log(s+1) - (k+1)*mp.log(s+2)
def dPhi(s,k):
    return j*pi/2 + (k-1)/s + 1/(s+1) - (k+1)/(s+2)
def ddPhi(s,k):
    return -(k-1)/s**2 - 1/(s+1)**2 + (k+1)/(s+2)**2

def Rmodel(k):
    rs=2*mp.sqrt(k/pi)
    def f(r):
        s=r*omega
        return H(s)*mp.expj(pi*s/2)*s**(k-1)*(s+1)/(s+2)**(k+1)*omega
    pts=[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf]
    return -4/pi*mp.quad(f,pts).imag

def saddle(k):
    r0=2*mp.sqrt((k-1)/pi)
    s0=r0*omega
    s_=mp.findroot(lambda s: dPhi(s,k), s0)
    GM = -(4/pi)*( omega*mp.exp(Phi(s_,k))*H(s_)*mp.sqrt(2*pi/(-ddPhi(s_,k))) ).imag
    return s_, GM

print(" k     R_k/(-(1/N)e^{-sq2piN})        GM/(-(1/N)e^{-sq2piN})       R/GM      phase-fit")
vals=[]
for k in [200,300,400,800,1200,1600,2400]:
    N=k+1
    R=Rmodel(k); s_,GM=saddle(k)
    th=sq2pi*mp.sqrt(N)
    sc = -(mp.e**(-th))/N
    # fitted "constant" assuming cos(th+phi): value/sc = C cos(th+phi)
    vals.append((k,R,GM,th,sc))
    print(f"{k:5d}  {mp.nstr(R/sc,8):>12s}   {mp.nstr(GM/sc,8):>12s}   {mp.nstr(R/GM,6)}   arg(Rat)={mp.nstr(mp.arg(R/sc),6)}")

# Now extract C, phi from GM: GM/sc should -> C cos(th+phi) with C,phi from saddle formula
# Instead print GM * N e^{th} as complex number before Im: i.e. gk = (4/pi)*omega e^{Phi}H sqrt(2pi/-Phi'')*N e^{th}
print()
print(" k     gk := -(4/pi) omega e^{Phi(s*)} H(s*) sqrt(2pi/(-Phi'')) * N e^{sq(2piN)}   (complex)")
for k in [400,1600,6400,25600]:
    N=k+1
    s_,GM=saddle(k)
    th=sq2pi*mp.sqrt(N)
    gk = -(4/pi)*omega*mp.exp(Phi(s_,k))*H(s_)*mp.sqrt(2*pi/(-ddPhi(s_,k)))*N*mp.e**th
    print(f"{k:6d}   {mp.nstr(gk,12)}")
