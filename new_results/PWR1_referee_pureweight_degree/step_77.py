import mpmath as mp, pickle
mp.mp.dps = 70
Ic={}
def I(m,j):
    if (m,j) not in Ic: Ic[(m,j)]=mp.quad(lambda ss: ss**m*mp.log1p(mp.e**(-ss))**j,[0,0.5,2,8,32,mp.inf])
    return Ic[(m,j)]
e={}
with open('e_rm_mine.pkl','rb') as f:
    for key,val in pickle.load(f).items():
        r,m=map(int,key.split(',')); e[(r,m)]=mp.mpf(val)
def Dt(n,r): return 2*mp.fsum(e[(r,n-1-j)]*I(n-1-j,j)/mp.factorial(j) for j in range(1,n-2*r))
def Theta(n,u): return mp.fsum(Dt(n,m)*(-4*u)**m/mp.factorial(2*m) for m in range(0,n//2))
def a_k(k,tau):
    tau=mp.mpf(tau); q=1/tau
    return mp.mpf(4)/mp.pi*mp.quad(lambda th:(mp.cos(th)**q+mp.sin(th)**q)**tau*mp.cos(4*k*th),[0,mp.pi/4,mp.pi/2])
taus=[mp.mpf(x) for x in ['0.001','0.002','0.003','0.004','0.006']]
Ak={k:-(-1)**k*4*mp.sqrt(2)/(mp.pi*(16*k*k-1)) for k in range(5)}
ak=[[a_k(k,t) for t in taus] for k in range(5)]
for n in range(2,6):
    mx=0
    for k in range(5):
        A=mp.matrix(len(taus),5); b=mp.matrix(len(taus),1)
        for i,tt in enumerate(taus):
            for j,nn in enumerate(range(2,7)): A[i,j]=tt**nn
            b[i]=ak[k][i]-Ak[k]
        sol=mp.lu_solve(A.T*A,A.T*b)
        v=(-1)**k*mp.pi/mp.sqrt(2)*sol[n-2]
        mx=max(mx,abs(v-Theta(n,k*k))/max(abs(Theta(n,k*k)),mp.mpf('1e-30')))
    print(f"fresh 5-pt extraction, n={n}: max rel.err over k=0..4 = {mp.nstr(mx,2)}")
