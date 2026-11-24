import time, sys, os, json
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
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

# Compute M values on circle ONCE, cache, then extract all k
def circle_vals(a,P,r,tag):
    fn=f'cache_{tag}.json'
    if os.path.exists(fn):
        d=json.load(open(fn)); return [mp.mpc(x[0],x[1]) for x in d]
    vals=[None]*P
    for p in range(P//2+1):
        th=2*mp.pi*p/P; u=r*mp.e**(1j*th)
        vals[p]=M_split2(u,a)
    for p in range(P//2+1,P): vals[p]=mp.conj(vals[P-p])
    json.dump([[float(mp.re(v)),float(mp.im(v))] for v in vals],open(fn,'w'))
    return vals

def extract(vals,Kmax,P,r):
    out=[]
    for k in range(Kmax+1):
        s=mp.fsum(vals[p]*mp.e**(-2j*mp.pi*k*p/P) for p in range(P))
        out.append(s/(P*r**k))
    return out

P=256; r=mp.mpf('0.9'); a=mp.pi/2
t=time.time()
vals=circle_vals(a,P,r,'pi2_P256')
cff=extract(vals,120,P,r)
print("time:",round(time.time()-t,1),"s")
# residual ratio at pi/2: |M_k + sqrt(pi)/N e^{-sqrt(pi N)} cos sqrt(pi N)| / (N^{-3/2} e^{-sqrt(pi N)})
def resid_ratio(Mk,k,a):
    N=k+1; A=mp.sqrt(2*a*mp.sin(a)); th=mp.sqrt(2*a*N)
    lead=-(A/N)*mp.e**(-th)*mp.cos(th)
    sc=N**(-mp.mpf(3)/2)*mp.e**(-th)
    return abs(Mk-lead)/sc, A, th
for k in [30,60,100,120]:
    rr,A,th=resid_ratio(mp.re(cff[k]),k,a)
    print(f"k={k}: M_k={mp.nstr(mp.re(cff[k]),10)} resid_ratio={mp.nstr(rr,6)}")
