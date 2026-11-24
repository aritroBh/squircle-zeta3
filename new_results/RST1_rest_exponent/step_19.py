"""RST1 step_19 (Task item 1): which singularity produces the rest + leading asymptotics.
(1) rest = M_k - M_main matches the q=2 model R_k = -(4/pi) Im int H E s^{k-1}(s+1)/(s+2)^{k+1} ds
    (the pole at s=-2 of the weight-kernel, i.e. the m=2 image of the +-i pi/2 branch points).
(2) Saddle of R_k: amplitude -> sqrt(pi)/N, phase -> theta_2 = sqrt(2 pi N).
(3) Subdominance of the q>=3 poles and the K-E (i.e. -E^3, s=3) image.
"""
import mpmath as mp, re
mp.mp.dps = 80
pi=mp.pi; j=1j; omega=mp.expjpi(mp.mpf(1)/4)
def H(s): return mp.beta(s/2, mp.mpf(3)/2)/2
sq2pi=mp.sqrt(2*pi); sqpi=mp.sqrt(pi)
def load():
    out={}
    for line in open("<project>/research/08_theory/D13f_sweep/mk_deg9.tsv"):
        p=line.rstrip("\n").split("\t")
        if not p[0].isdigit(): continue
        out[int(p[0])]=mp.mpf(re.match(r"\[([-0-9.eE+]+)",p[1]).group(1))
    return out
rows=load()
def Mmain(k):
    N=k+1; rs=mp.sqrt(2*N/pi)
    f=lambda r: mp.expj(pi*(r*omega)/2)*(r*omega/(1+r*omega))**N*H(r*omega)/(r*omega)*omega
    return 2/pi*mp.quad(f,[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf]).imag
def Rq(k,q):
    """model integral for pole at s=-q: c_q Im int H E s^{k-q+1}(s+1)^{q-1}/(s+q)^{k+1} ds (sign/structure from A_k expansion)"""
    rs=mp.sqrt(2*k*q/pi)
    def f(r):
        s=r*omega
        return H(s)*mp.expj(pi*s/2)*s**(k-1)*(s+1)/(s+2)**(k+1)*omega if q==2 else None
    return None
def Rmodel2(k):
    rs=2*mp.sqrt(k/pi)
    def f(r):
        s=r*omega
        return H(s)*mp.expj(pi*s/2)*s**(k-1)*(s+1)/(s+2)**(k+1)*omega
    return -4/pi*mp.quad(f,[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf]).imag
def Rmodel3(k):
    """q=3 pole model: leading term (1/2) s^{k-2}/(s+3)^{k-1} (from i=2 of A_k), sign +"""
    rs=mp.sqrt(6*k/pi)
    def f(r):
        s=r*omega
        return H(s)*mp.expj(pi*s/2)*s**(k-2)/(s+3)**(k-1)*omega
    return (1/pi)*mp.quad(f,[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf]).imag
def RKE3(k):
    """K-E image: -E^3 Q_k model: -(2/pi) Im int H (-E^3) s^k/(1+s)^{k+1} ds"""
    rs=mp.sqrt(2*k/(3*pi))
    def f(r):
        s=r*omega
        return H(s)*mp.expj(3*pi*s/2)*s**k/(1+s)**(k+1)*omega
    return (2/pi)*mp.quad(f,[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf]).imag

print("== (1) rest vs q=2 model R_k (exact M_k) ==")
print(" k      rest               R_k              rest/R_k")
for k in [100,150,200,250,300,350,400]:
    rest=rows[k]-Mmain(k); Rm=Rmodel2(k)
    print(f"{k:4d}  {mp.nstr(rest,8):>14s}  {mp.nstr(Rm,8):>14s}   {mp.nstr(rest/Rm,7)}")

print()
print("== (2) saddle leading law of R_k:  R_k ~ (sqrt pi/N) e^{-theta_2} cos(theta_2) ==")
print(" k         R_k/lead          |gk|*N e^{th2} (->sqrt pi)   arg(gkf)-th2 (->pi/2)")
def saddle(k):
    r0=2*mp.sqrt((k-1)/pi)
    dPhi=lambda s: j*pi/2+(k-1)/s+1/(s+1)-(k+1)/(s+2)
    Phi=lambda s: j*pi*s/2+(k-1)*mp.log(s)+mp.log(s+1)-(k+1)*mp.log(s+2)
    ddPhi=lambda s: -(k-1)/s**2-1/(s+1)**2+(k+1)/(s+2)**2
    s_=mp.findroot(dPhi, r0*omega)
    g=omega*mp.exp(Phi(s_))*H(s_)*mp.sqrt(2*pi/(-ddPhi(s_)))
    return -(4/pi)*g
for k in [400,2000,10000,50000]:
    N=k+1; th=sq2pi*mp.sqrt(N)
    gkf=saddle(k); Rm=Rmodel2(k)
    lead=sqpi/N*mp.e**(-th)*mp.cos(th)
    print(f"{k:6d}  {mp.nstr(Rm/lead,8):>10s}    {mp.nstr(abs(gkf)*N*mp.e**th,10):>12s}        {mp.nstr(mp.fmod(mp.arg(gkf)-th+pi,2*pi)-pi,10)}")

print()
print("== (3) subdominance: |q=3 model| and |K-E model| vs |q=2 model| ==")
print(" k      |R_q3|/|R_q2|      |R_KE3|/|R_q2|")
for k in [100,200,400]:
    r2=abs(Rmodel2(k)); r3=abs(Rmodel3(k)); rke=abs(RKE3(k))
    print(f"{k:4d}    {mp.nstr(r3/r2,6):>10s}      {mp.nstr(rke/r2,6)}")
