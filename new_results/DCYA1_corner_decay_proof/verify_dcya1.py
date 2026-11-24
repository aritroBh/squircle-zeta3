"""verify_dcya1.py — checks the PROVED exact representation (L1) and the main-term structure (L4)
for DCYA1. Run: <venv>/python verify_dcya1.py   (should finish in a few minutes)"""
import mpmath as mp
mp.mp.dps=50
def gm(c,M):
    Q=[mp.mpf(0)]*(M+1); Q[0]=1
    if M>=1:Q[1]=-2*c
    if M>=2:Q[2]=1
    b=[mp.sqrt(Q[0])]
    for n in range(1,M+1):
        b.append((Q[n]-mp.fsum(b[k]*b[n-k] for k in range(1,n)))/(2*b[0]))
    return b
def bj(u,J):
    o=[mp.mpf(1)]
    for j in range(1,J+1): o.append(-(u+j)/j*o[-1])
    return o
def M_def(u,a,Jm=20):
    c=mp.cos(a); g0=mp.mpf(1)/64; G=gm(c,Jm); B=bj(u,Jm); s=mp.mpf(0)
    for j in range(Jm+1):
        for m in range(Jm+1):
            if(j,m)!=(0,0):
                nu=j+(m-1)*(1-u); s+=B[j]*G[m]*g0**(nu+1)/(nu+1)
    F=lambda g:((1+g)**(-1-u)*mp.sqrt(1+g**(2-2*u)-2*c*g**(1-u))-1)*g**(u-1)
    return s+mp.quad(F,[g0,mp.mpf(1)/8,mp.mpf(1)/2,1])
def Bt(z,a):
    c=mp.cos(a); return (mp.e**z-2*c)/(1+mp.sqrt(1-2*c*mp.e**z+mp.e**(2*z)))
def T_u(u,a):
    L=lambda x:mp.log(1+mp.e**(-x)); return mp.quad(lambda x:mp.e**(-(1+u)*L(x))*mp.e**(-x)*Bt((u-1)*x,a),[0,1,4,16,mp.inf])
def R_u(u,a):
    L=lambda x:mp.log(1+mp.e**(-x)); return mp.quad(lambda x:(mp.e**(-(1+u)*L(x))-1)*mp.e**(-u*x),[0,1,4,16,mp.inf])
print("L1: |M(u;a)-(T+R)| (exact representation)")
for a in [mp.pi/3,2*mp.pi/3,3*mp.pi/4,mp.pi/2]:
    for u in [mp.mpf('0.1'),mp.mpf('0.3')]:
        print(f"  a={mp.nstr(a,4)} u={u}: {mp.nstr(abs(M_def(u,a)-(T_u(u,a)+R_u(u,a))),3)}")
# L4 structure check: main term leading+g_1 captures kappa_1 at small N via direct I_+ integrals
def Iplus(k,a,jj):
    N=k+1; E=mp.e**(1j*a); gp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))
    g1=-mp.mpf(3)/4-1j*E/(4*mp.sin(a)); gj=[1,g1][jj]
    xs=mp.sqrt(N*a)
    J=mp.quad(lambda x:mp.e**(-x)*x**(N-1)*(x+1j*a)**(-(N-jj-mp.mpf(3)/2)),[0,xs/4,xs,4*xs,mp.inf])
    return gp*gj*(1j/mp.pi)*mp.gamma(jj+mp.mpf(3)/2)*mp.gamma(N-jj-mp.mpf(3)/2)/mp.gamma(N)*J
print("L4: main term leading+g_1 ~ M_k (ratio ~1+O(1/sqrt N)); uses quad M_k via series-evaluator")
import corner_series as cs
for a in [2*mp.pi/3]:
    k=40
    Mk=cs.corner_coeffs(a,k)[k]
    mm=2*mp.re(Iplus(k,a,0)+Iplus(k,a,1))
    print(f"  a=2pi/3 k={k}: M_k={mp.nstr(Mk,10)}  main(leading+g_1)={mp.nstr(mm,10)}  ratio={mp.nstr(Mk/mm,5)}")
print("done.")
