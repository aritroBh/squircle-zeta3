import time, sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
import corner_series as cs
def gamma_m(c, M):
    Q=[mp.mpf(0)]*(M+1); Q[0]=1
    if M>=1: Q[1]=-2*c
    if M>=2: Q[2]=1
    b=[mp.sqrt(Q[0])]
    for n in range(1,M+1):
        s=Q[n]-mp.fsum(b[k]*b[n-k] for k in range(1,n))
        b.append(s/(2*b[0]))
    return b
def beta_j(u, J):
    out=[mp.mpf(1)]
    for j in range(1,J+1): out.append(-(u+j)/j*out[-1])
    return out
def M_split2(u, a, Jm=30, g0=None):
    c=mp.cos(a)
    if g0 is None: g0=mp.mpf(1)/64
    gm=gamma_m(c,Jm); bj=beta_j(u,Jm)
    ana=mp.mpf(0)
    for j in range(Jm+1):
        for m in range(Jm+1):
            if (j,m)==(0,0): continue
            nu=j+(m-1)*(1-u)
            ana+=bj[j]*gm[m]*g0**(nu+1)/(nu+1)
    def F(g):
        return ((1+g)**(-1-u)*mp.sqrt(1+g**(2-2*u)-2*c*g**(1-u))-1)*g**(u-1)
    return ana+mp.quad(F, [g0, mp.mpf(1)/8, mp.mpf(1)/2, 1])

def Mk_fft(a, Kmax, P=256, r=mp.mpf('0.9')):
    vals=[None]*P
    for p in range(P):
        th=2*mp.pi*p/P
        u=r*mp.e**(1j*th)
        if p<=P//2: vals[p]=M_split2(u,a)
        else: vals[p]=mp.conj(vals[P-p])
    out=[]
    for k in range(Kmax+1):
        s=mp.fsum(vals[p]*mp.e**(-2j*mp.pi*k*p/P) for p in range(P))
        out.append(s/(P*r**k))
    return out

t=time.time()
a=mp.pi/2
K=25
cff=Mk_fft(a,K,P=256)
print("FFT time:",round(time.time()-t,1),"s")
cref=cs.corner_coeffs(a,20)
for k in [0,1,2,5,10,15,20]:
    err=abs(cff[k]-cref[k])
    print(f"k={k}: |err|={mp.nstr(err,3)}  val={mp.nstr(mp.re(cff[k]),15)}")
