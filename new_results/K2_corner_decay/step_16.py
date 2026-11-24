
import mpmath as mp, time
mp.mp.dps = 100
P=1024; KMAX=300; R=mp.mpf('0.9'); M1=J1=6

def gamma_coeffs(c, M1):
    Q=[mp.mpf(1),-2*c,mp.mpf(1)]+[mp.mpf(0)]*(M1-2)
    b=[mp.sqrt(Q[0])]
    for n in range(1,M1+1):
        s=Q[n]-mp.fsum(b[k]*b[n-k] for k in range(1,n)); b.append(s/(2*b[0]))
    return b

def make_grid(H,JTS):
    grid=[]
    for jj in range(-JTS,JTS+1):
        t=jj*H; st=mp.sinh(t); th=mp.tanh(mp.pi/2*st); g=(1+th)/2
        w=(mp.pi/4)*mp.cosh(t)/(mp.cosh(mp.pi/2*st)**2)
        if g<=0 or g>=1: continue
        grid.append((g,w,mp.log(g),mp.log1p(g)))
    return grid

alpha=mp.pi/2; c=mp.cos(alpha); gam=gamma_coeffs(c,M1)
def Mval(u, grid, H):
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

grid_c = make_grid(mp.mpf('0.05'), 112)   # coarse h=0.05
grid_f = make_grid(mp.mpf('0.025'), 224)  # fine h=0.025
t0=time.time()
errs=[]
for p in range(0, P//2+1):
    u=R*mp.exp(2*mp.pi*mp.j*p/P)
    vc=Mval(u, grid_c, mp.mpf('0.05'))
    vf=Mval(u, grid_f, mp.mpf('0.025'))
    errs.append((abs(vf-vc), p))
errs.sort(reverse=True)
print("scan done %.1fs; 10 worst (|fine-coarse|, p):" % (time.time()-t0))
for e,p in errs[:10]:
    u=R*mp.exp(2*mp.pi*mp.j*p/P)
    print("  p=%3d u=%.3f%+.3fi  diff=%.2e" % (p, u.real, u.imag, e))
