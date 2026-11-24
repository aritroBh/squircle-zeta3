import time, sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps = 60

def gamma_m(c, M):
    """coeffs [h^m] sqrt(1-2 c h + h^2), m=0..M"""
    # series sqrt
    Q=[mp.mpf(0)]*(M+1); Q[0]=1
    if M>=1: Q[1]=-2*c
    if M>=2: Q[2]=1
    b=[mp.sqrt(Q[0])]
    for n in range(1,M+1):
        s=Q[n]-mp.fsum(b[k]*b[n-k] for k in range(1,n))
        b.append(s/(2*b[0]))
    return b

def beta_j(u, J):
    """coeffs [g^j] (1+g)^(-1-u) = (-1)^j (1+u)_j / j!"""
    out=[mp.mpf(1)]
    for j in range(1,J+1):
        out.append(-(u+j)/j*out[-1])
    return out

# Verify the small-g expansion against the actual integrand
c=mp.cos(mp.pi/2)
u=mp.mpf('-0.9')
Jm=12
gm=gamma_m(c,Jm); bj=beta_j(u,Jm)
g=mp.mpf('0.01')
# full integrand
Ffull=((1+g)**(-1-u)*mp.sqrt(1+g**(2-2*u)-2*c*g**(1-u))-1)*g**(u-1)
# series
Fser=mp.fsum(bj[j]*gm[m]*g**(j+(m-1)*(1-u)) for j in range(Jm+1) for m in range(Jm+1) if (j,m)!=(0,0))
print("Ffull=",mp.nstr(Ffull,20))
print("Fser =",mp.nstr(Fser,20))
print("diff=",mp.nstr(abs(Ffull-Fser),4))
print("gamma_1..5:",[mp.nstr(x,8) for x in gm[:6]])
