import time, sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
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
    for j in range(1,J+1):
        out.append(-(u+j)/j*out[-1])
    return out

def M_sub(u, a, Jm=14, pts=None):
    c=mp.cos(a)
    gm=gamma_m(c,Jm); bj=beta_j(u,Jm)
    terms=[(j,m) for j in range(Jm+1) for m in range(Jm+1) if (j,m)!=(0,0)]
    def Fsub(g):
        Ffull=((1+g)**(-1-u)*mp.sqrt(1+g**(2-2*u)-2*c*g**(1-u))-1)*g**(u-1)
        S=mp.fsum(bj[j]*gm[m]*g**(j+(m-1)*(1-u)) for (j,m) in terms)
        return Ffull-S
    if pts is None:
        pts=[0, mp.mpf(1)/64, mp.mpf(1)/8, mp.mpf(1)/2, 1]
    I=mp.quad(Fsub, pts)
    analytic=mp.fsum(bj[j]*gm[m]/(j+(m-1)*(1-u)+1) for (j,m) in terms)
    return I+analytic

# validate at u=-0.9, pi/2 against high-precision direct
mp.mp.dps=60
for dps in [40,60,80]:
    mp.mp.dps=dps
    v=M_sub(mp.mpf('-0.9'), mp.pi/2, Jm=14)
    print(f"dps={dps}: M_sub(-0.9)={mp.nstr(v,32)}")
mp.mp.dps=120
vref=cs.M_exact(mp.mpf('-0.9'), mp.pi/2)
print("ref(dps120,direct)=",mp.nstr(vref,32))
