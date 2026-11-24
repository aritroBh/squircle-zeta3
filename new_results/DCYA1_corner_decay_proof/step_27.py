import time, sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
import corner_series as cs
exec(open('/dev/stdin').read()) if False else None
# redefine with tunable quad degree
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
def M_split2(u, a, Jm=18, g0=None, deg=None):
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
    old=mp.mp.prec
    I=mp.quad(F, [g0, mp.mpf(1)/8, mp.mpf(1)/2, 1], degree=deg) if deg else mp.quad(F, [g0, mp.mpf(1)/8, mp.mpf(1)/2, 1])
    return ana+I
mp.mp.dps=100
vref=cs.M_exact(mp.mpf('-0.9'), mp.pi/2)
mp.mp.dps=60
for Jm in [18,30,50]:
    t=time.time(); v=M_split2(mp.mpf('-0.9'),mp.pi/2,Jm=Jm); 
    print(f"Jm={Jm}: err={mp.nstr(abs(v-vref),3)} t={round(time.time()-t,2)}s")
# try higher quad degree
v=M_split2(mp.mpf('-0.9'),mp.pi/2,Jm=30,deg=6)
print("Jm=30 deg=6 err=",mp.nstr(abs(v-vref),3))
