
# step_10.py : fix I44, I26 relations (tight PSLQ); independent d9 check via series Richardson;
# consolidated square-end table.
import mpmath as mp
mp.mp.dps = 130
import sympy as sp
Kmax = 70
u = sp.symbols('u')
ex = sp.log(sum((-u)**j/sp.factorial(j+1) for j in range(Kmax+3)))
cs=[sp.expand(ex).series(u,0,Kmax+1).removeO().coeff(u,k) for k in range(Kmax+1)]
c=[mp.mpf(x) for x in cs]
def gpow(j):
    poly=[mp.mpf(1)]
    for _ in range(j):
        new=[mp.mpf(0)]*(len(poly)+Kmax)
        for i,a in enumerate(poly):
            for k2 in range(1,Kmax+1): new[i+k2]+=a*c[k2]
        poly=new
    return poly
GP={j:gpow(j) for j in range(0,11)}
def int_u_log(s,m,delta):
    ld=mp.log(delta); tot=mp.mpf(0)
    for i in range(m+1): tot+=(-1)**(m-i)*mp.factorial(m)/mp.factorial(i)*ld**i/(s+1)**(m-i+1)
    return delta**(s+1)*tot
def I_ab(a,b,delta=mp.mpf('0.5')):
    near=mp.mpf(0)
    for j in range(0,b+1):
        for k2,coeff in enumerate(GP[j]):
            if coeff==0: continue
            near+=mp.binomial(b,j)*coeff*int_u_log(a+k2,b-j,delta)
    return near+mp.quad(lambda uu: uu**a*mp.log(1-mp.exp(-uu))**b,[delta,1,3,8,20,mp.inf])
zb=[mp.zeta(9),mp.zeta(2)*mp.zeta(7),mp.pi**6*mp.zeta(3),mp.pi**4*mp.zeta(5),mp.zeta(3)**3]
names=['z9','z2z7','pi^6z3','pi^4z5','z3^3']
for ab in [(4,4),(2,6)]:
    v=I_ab(*ab)
    rel=mp.pslq([v]+zb,tol=mp.mpf('1e-70'),maxcoeff=10**13,maxsteps=10000)
    print(f'I_{ab} tight PSLQ:',rel)
    if rel:
        q=rel[0]; val=sum(-mp.mpf(rel[i+1])/q*b for i,b in enumerate(zb))
        print('   residual:',mp.nstr(abs(v-val),4),' coeffs:',[(str(sp.Rational(rel[i+1],q)),names[i]) for i in range(5)])

# independent d9 via series Richardson extraction (10-digit confirmation)
mp.mp.dps = 105
def A_series(p,tol=mp.mpf('1e-101')):
    p=mp.mpf(p); s=mp.mpf(0); m=0
    while True:
        m+=1
        t=mp.gamma(1+m/p)**2/(m**2*mp.gamma(1+2*m/p))
        s+=t
        if m>20 and abs(t)<tol: break
    return s
ps=list(range(60,501,10))
xs=[mp.mpf(1)/p for p in ps]
ys=[mp.zeta(2)-A_series(p) for p in ps]
deg=10
n=len(xs); M=mp.matrix(n,deg)
for i,x in enumerate(xs):
    for j in range(deg): M[i,j]=x**(j+1)
b=mp.matrix(n,1)
for i,y in enumerate(ys): b[i,0]=y
sol=mp.lu_solve(M.T*M,M.T*b)
print('series-extracted d1..d10:',[mp.nstr(sol[j,0],11) for j in range(deg)])
d9exact=(-428400*zb[0]+136080*zb[1]+79*zb[2]+1134*zb[3]-10080*zb[4])/1814400
print('d9 exact candidate =', mp.nstr(d9exact,25))
