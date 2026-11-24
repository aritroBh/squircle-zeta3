
import mpmath as mp, time, json
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

Hf=mp.mpf('0.025'); grid_f=make_grid(Hf,224)
t0=time.time()
f=[Mval(R*mp.exp(2*mp.pi*mp.j*p/P), grid_f, Hf) for p in range(P//2+1)]
print("f computed %.1fs"%(time.time()-t0))
json.dump([[str(x.real),str(x.imag)] for x in f], open('out/fval_pi2_h025.json','w'))

half=P//2; om=mp.exp(-2*mp.pi*mp.j/P)
ims=[]
for k in range(KMAX+1):
    s=mp.mpc(0); wk=om**k; wp=mp.mpc(1)
    for p in range(1,half):
        wp*=wk; s+=f[p]*wp
    ims.append((abs(2*s.imag/P), k))
ims.sort(reverse=True)
print("worst |Im a_k|:", [("%.2e"%e, k) for e,k in ims[:6]])
print("f max abs:", max(abs(x) for x in f))
print("f[0]=%s\nf[512]=%s"%(mp.nstr(f[0],20), mp.nstr(f[512],20)))
