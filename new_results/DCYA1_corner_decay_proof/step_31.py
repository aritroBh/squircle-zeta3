import time, sys, os
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=60
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
vals=[None]*P
for p in range(P//2+1):
    th=2*mp.pi*p/P; u=r*mp.e**(1j*th); vals[p]=M_split2(u,a)
for p in range(P//2+1,P): vals[p]=mp.conj(vals[P-p])
# DFT sum for k=100: magnitude before dividing
for k in [20,60,100]:
    s=mp.fsum(vals[p]*mp.e**(-2j*mp.pi*k*p/P) for p in range(P))
    Mk=s/(P*r**k)
    print(f"k={k}: |sum|={mp.nstr(abs(s),4)}  M_k={mp.nstr(mp.re(Mk),8)}  max|val|={mp.nstr(max(abs(v) for v in vals),4)}")
# Check: what is the quadrature noise level? Estimate via imaginary part of M_k (should be ~0 for real M_k)
print("M_k are real; check imaginary parts:")
for k in [20,60,100]:
    s=mp.fsum(vals[p]*mp.e**(-2j*mp.pi*k*p/P) for p in range(P))
    Mk=s/(P*r**k)
    print(f"k={k}: im(M_k)={mp.nstr(mp.im(Mk),3)}")
