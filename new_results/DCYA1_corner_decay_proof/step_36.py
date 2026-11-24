import time, sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
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
def M_split2(u, a, Jm=30):
    c=mp.cos(a); g0=mp.mpf(1)/64
    gm=gamma_m(c,Jm); bj=beta_j(u,Jm)
    ana=mp.mpf(0)
    for j in range(Jm+1):
        for m in range(Jm+1):
            if (j,m)==(0,0): continue
            nu=j+(m-1)*(1-u); ana+=bj[j]*gm[m]*g0**(nu+1)/(nu+1)
    def F(g):
        return ((1+g)**(-1-u)*mp.sqrt(1+g**(2-2*u)-2*c*g**(1-u))-1)*g**(u-1)
    return ana+mp.quad(F, [g0, mp.mpf(1)/8, mp.mpf(1)/2, 1])

P=256; r=mp.mpf('0.9'); a=mp.pi/2
for dps in [50,90]:
    mp.mp.dps=dps
    vals=[None]*P
    for p in range(P//2+1):
        u=r*mp.e**(2j*mp.pi*p/P); vals[p]=M_split2(u,a)
    for p in range(P//2+1,P): vals[p]=mp.conj(vals[P-p])
    for k in [100]:
        s=mp.fsum(vals[p]*mp.e**(-2j*mp.pi*k*p/P) for p in range(P))
        Mk=s/(P*r**k)
        print(f"dps={dps} k={k}: M_k={mp.nstr(mp.re(Mk),16)}  |sum|={mp.nstr(abs(s),4)}")
