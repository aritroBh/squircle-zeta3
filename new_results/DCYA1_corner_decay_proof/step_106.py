import sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=40
b=mp.e**(1j*mp.pi/4)
def Phi(t,N,a,j=0): return -b*t+(N-1)*mp.log(b*t)-(N-j-mp.mpf(3)/2)*mp.log(b*t+1j*a)
def Phip(t,N,a,j=0): return -b+(N-1)/t-(N-j-mp.mpf(3)/2)*b/(b*t+1j*a)
def Phipp(t,N,a,j=0): return -(N-1)/t**2+(N-j-mp.mpf(3)/2)*b**2/(b*t+1j*a)**2
def gmoment0(p2,om): return mp.sqrt(2*mp.pi/(-p2))*mp.e**(om**2/(2*p2))
a=2*mp.pi/3; A=mp.sqrt(2*a*mp.sin(a)); N=400
t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,0)), mp.sqrt(N*a))
om=mp.im(Phip(t0,N,a,0)); p2=Phipp(t0,N,a,0); ps=Phi(t0,N,a,0)
pre=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)
MV0=pre*mp.e**ps*gmoment0(p2,om)*b
th=mp.sqrt(2*a*N); sc=N**(-1.5)*mp.e**(-th); lead=-(A/N)*mp.e**(-th)*mp.cos(th)
print("mpmath: 2ReMV0=",mp.nstr(2*mp.re(MV0),8)," lead=",mp.nstr(lead,8))
print("mpmath: |2ReMV0-lead|/sc=",mp.nstr(abs(2*mp.re(MV0)-lead)/sc,4)," sc=",mp.nstr(sc,4))
print("cos(th)=",mp.nstr(mp.cos(th),6))
