
import mpmath as mp
mp.mp.dps = 100

# reload f-values? They weren't saved. Recompute grid + f for pi/2 quickly at dps 100, then diagnose.
import time
P=1024; KMAX=300; R=mp.mpf('0.9'); H=mp.mpf('0.025'); JTS=224; M1=J1=6

def gamma_coeffs(c, M1):
    Q=[mp.mpf(1),-2*c,mp.mpf(1)]+[mp.mpf(0)]*(M1-2)
    b=[mp.sqrt(Q[0])]
    for n in range(1,M1+1):
        s=Q[n]-mp.fsum(b[k]*b[n-k] for k in range(1,n)); b.append(s/(2*b[0]))
    return b

grid=[]
for jj in range(-JTS,JTS+1):
    t=jj*H; st=mp.sinh(t); th=mp.tanh(mp.pi/2*st); g=(1+th)/2
    w=(mp.pi/4)*mp.cosh(t)/(mp.cosh(mp.pi/2*st)**2)
    if g<=0 or g>=1: continue
    grid.append((g,w,mp.log(g),mp.log1p(g)))

alpha=mp.pi/2; c=mp.cos(alpha); gam=gamma_coeffs(c,M1)
def Mval(u):
    beta=[mp.mpc(1)]
    for j in range(J1): beta.append(beta[-1]*(-(1+u+j))/(j+1))
    ana=mp.mpc(0)
    for m in range(1,M1+1): ana+=gam[m]/(m-(m-1)*u)
    for j in range(1,J1+1): ana+=beta[j]/(u+j)
    tot=mp.mpc(0)
    for (g,w,lg,l1g) in grid:
        A=mp.exp(lg*(u-1)); B=1/A; lC=-(1+u)*l1g; C=mp.exp(lC)
        Q=1-2*c*B+B*B; sQ=mp.sqrt(Q)
        Cm1=C*((Q-1)/(sQ+1))+mp.expm1(lC)
        S=Cm1*A
        Bp=mp.mpc(1)
        for m in range(1,M1+1): S-=gam[m]*Bp; Bp*=B
        Ap=A*g
        for j in range(1,J1+1): S-=beta[j]*Ap; Ap*=g
        tot+=w*S
    return tot*H+ana

# direct high-accuracy reference via mp.quad on original integrand
def Mref(u):
    f=lambda g: ((1+g)**(-1-u)*mp.sqrt(1+g**(2-2*u)-2*c*g**(1-u))-1)*g**(u-1)
    return mp.quad(f,[0,mp.mpf('0.01'),mp.mpf('0.5'),mp.mpf('0.99'),1])

worst=[]
for p in [0, 64, 128, 192, 256, 320, 384, 448, 500, 512]:
    u=R*mp.exp(2*mp.pi*mp.j*p/P)
    v1=Mval(u)
    if p<=300:
        v2=Mref(u)
        err=abs(v1-v2)
    else:
        err=None
    print("p=%3d u=%.3f%+.3fi  |Mval|=%.3e  err_vs_quad=%s" % (p, u.real, u.imag, abs(v1), ("%.2e"%err) if err else "skip"))
