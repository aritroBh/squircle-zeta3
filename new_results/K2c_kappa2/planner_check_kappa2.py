# Planner check (run from research/08_theory/K2b_corner_phase): kappa2 formula vs corner M_k data, k <= 300, five angles.
import mpmath as mp
mp.mp.dps=40
angles={'pi_3':mp.pi/3,'pi_2':mp.pi/2,'2pi_3':2*mp.pi/3,'3pi_4':3*mp.pi/4,'5pi_6':5*mp.pi/6}
def load(tag):
    M={}
    for line in open(f'../K2_corner_decay/out/mk_{tag}.tsv'):
        if line.startswith('#') or not line.strip(): continue
        k,v=line.split('\t'); t=v.strip().strip('()').replace(' ','')[:-1]
        i=max(j for j in range(1,len(t)) if t[j] in '+-' and t[j-1] not in 'eE'); M[int(k)]=mp.mpf(t[:i])
    return M
def kappa1(a):
    C=1/(2*mp.sqrt(a))+a**1.5/12+3*mp.sqrt(a)/8*mp.cot(a)
    return (C+mp.sqrt(a)/2)/mp.sqrt(2)+1j*(C-mp.sqrt(a)/2)/mp.sqrt(2)
def kappa2(a):
    ct=mp.cot(a); cs2=1/mp.sin(a)**2
    return (mp.mpf(3)/4+a**2/24+3*a/16*ct) + 1j*(a**3/288-11*a/384+(a**2/32+mp.mpf(9)/16)*ct-15*a/128*cs2)
print("kappa2(pi/2) =", mp.nstr(kappa2(mp.pi/2),20), " closed:", mp.nstr(mp.mpf(3)/4+mp.pi**2/96 - 1j*(7*mp.pi/96-mp.pi**3/2304),20))
for tag,a in angles.items():
    M=load(tag); A=mp.sqrt(2*a*mp.sin(a)); k1=kappa1(a); k2=kappa2(a)
    # least-squares estimate of kappa2 from data: N*rho1 ~ -(cos th Re k2 - sin th Im k2)
    rows=[];rhs=[]
    res={}
    for lab,kk in [('k1 only',0),('k1+k2',k2),('k1+conj k2',mp.conj(k2))]:
        mx=0
        for k in range(150,max(M)+1):
            N=k+1; th=mp.sqrt(2*a*N)
            pred=-(A/N)*mp.e**(-th)*mp.re(mp.e**(1j*th)*(1+k1/mp.sqrt(N)+kk/N))
            rho=(M[k]-pred)*N*mp.e**th/A
            mx=max(mx,abs(rho)*(N**1.5 if kk!=0 else N))
        res[lab]=mx
    for k in range(200,max(M)+1):
        N=k+1; th=mp.sqrt(2*a*N)
        pred=-(A/N)*mp.e**(-th)*mp.re(mp.e**(1j*th)*(1+k1/mp.sqrt(N)))
        y=(M[k]-pred)*N*mp.e**th/A*N
        rows.append([-mp.cos(th), mp.sin(th), -mp.cos(th)/mp.sqrt(N), mp.sin(th)/mp.sqrt(N)]); rhs.append(y)
    X=mp.matrix(rows); Y=mp.matrix(rhs)
    sol=mp.lu_solve(X.T*X, X.T*Y)
    print(f"{tag:6s} max|rho|*N (k1 only)={mp.nstr(res['k1 only'],4)}  max|rho|*N^1.5: k2={mp.nstr(res['k1+k2'],4)}, conj={mp.nstr(res['k1+conj k2'],4)} | fit k2 = {mp.nstr(sol[0],8)}{'+' if sol[1]>=0 else ''}{mp.nstr(sol[1],8)}i vs formula {mp.nstr(k2,8)}")
print("windowed max|rho|*N^1.5 with kappa2 (stable => kappa2 right; growth ~sqrt(N) => wrong); and with kappa2 perturbed by 0.05")
for tag,a in angles.items():
    M=load(tag); A=mp.sqrt(2*a*mp.sin(a)); k1=kappa1(a); k2=kappa2(a)
    def win(kk,lo,hi):
        mx=0
        for k in range(lo,hi):
            N=k+1; th=mp.sqrt(2*a*N)
            pred=-(A/N)*mp.e**(-th)*mp.re(mp.e**(1j*th)*(1+k1/mp.sqrt(N)+kk/N))
            mx=max(mx,abs((M[k]-pred)*N*mp.e**th/A)*N**1.5)
        return mx
    K=max(M)
    w=[(lo,min(lo+50,K+1)) for lo in (50,100,150,200,250)]
    print(f"{tag:6s} formula:", " ".join(mp.nstr(win(k2,lo,hi),3) for lo,hi in w), "| +0.05:", " ".join(mp.nstr(win(k2+0.05,lo,hi),3) for lo,hi in w), "| +0.05i:", " ".join(mp.nstr(win(k2+0.05j,lo,hi),3) for lo,hi in w))
