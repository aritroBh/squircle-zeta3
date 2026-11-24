import time, sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=70
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
def M_split4(u, a, Jm, g0):
    c=mp.cos(a)
    gm=gamma_m(c,Jm); bj=beta_j(u,Jm)
    ana=mp.mpf(0)
    for j in range(Jm+1):
        for m in range(Jm+1):
            if (j,m)==(0,0): continue
            nu=j+(m-1)*(1-u); ana+=bj[j]*gm[m]*g0**(nu+1)/(nu+1)
    def F(g):
        return ((1+g)**(-1-u)*mp.sqrt(1+g**(2-2*u)-2*c*g**(1-u))-1)*g**(u-1)
    return ana+mp.quad(F, [g0, (g0+1)/2, 1])
u=mp.mpf('-0.9'); a=mp.pi/2
for g0 in [mp.mpf(1)/16, mp.mpf(1)/64, mp.mpf(1)/256]:
    # need Jm large enough for small g0? use Jm=30
    v=M_split4(u,a,30,g0)
    print(f"g0=1/{int(1/g0)}: {mp.nstr(v,40)}")
