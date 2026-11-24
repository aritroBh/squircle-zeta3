import time, sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
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
def M_def(u, a, Jm=30):
    # M(u;a) via small-g split (analytic + numerical)
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

def Btilde(z, a):
    c=mp.cos(a)
    return (mp.e**z-2*c)/(1+mp.sqrt(1-2*c*mp.e**z+mp.e**(2*z)))
def T_of_u(u, a):
    L=lambda x: mp.log(1+mp.e**(-x))
    f=lambda x: mp.e**(-(1+u)*L(x))*mp.e**(-x)*Btilde((u-1)*x, a)
    return mp.quad(f,[0,1,4,16,mp.inf])
def R_of_u(u, a):
    L=lambda x: mp.log(1+mp.e**(-x))
    f=lambda x: (mp.e**(-(1+u)*L(x))-1)*mp.e**(-u*x)
    return mp.quad(f,[0,1,4,16,mp.inf])

print("Verify M(u;a) = T(u) + R(u):")
for a in [mp.pi/3, 2*mp.pi/3, 3*mp.pi/4, mp.pi/2]:
    for u in [mp.mpf('0.1'), mp.mpf('0.3')]:
        Md=M_def(u,a); T=T_of_u(u,a); R=R_of_u(u,a)
        err=abs(Md-(T+R))
        print(f"  a={mp.nstr(a,4)} u={u}: M={mp.nstr(Md,16)}  |M-(T+R)|={mp.nstr(err,3)}")
