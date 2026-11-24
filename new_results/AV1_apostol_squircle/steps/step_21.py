
# step_21.py : (i) d1 relation residual + 90dps rerun; (ii) d11 vs independent weight-11 basis; (iii) d11 Richardson.
import mpmath as mp
import sympy as sp

d1str = '1.73300259020287278417012270424573195250974925287880319940547987633644278103184087570877319'

for dps in [90, 130]:
    mp.mp.dps = dps
    d1 = mp.mpf(d1str)
    pool = [mp.mpf(1),mp.euler,mp.log(2*mp.pi),mp.log(2),mp.zeta(2),mp.zeta(3)]
    rel = mp.pslq([d1]+pool, tol=mp.mpf(10)**(-(dps-40)), maxcoeff=10**10, maxsteps=10000)
    print(f'dps={dps}: d1 pool1 PSLQ:', rel)
mp.mp.dps = 130
d1 = mp.mpf(d1str)
rel = [-2968836047, 446072542, 1981735866, 298434566, -878652262, -803572980, 4107470189]
pool = [mp.mpf(1),mp.euler,mp.log(2*mp.pi),mp.log(2),mp.zeta(2),mp.zeta(3)]
res = rel[0]*d1 + sum(r*b for r,b in zip(rel[1:],pool))
print('explicit residual of that relation:', mp.nstr(res, 6))

# d11 machinery
Kmax = 80
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
GP={j:gpow(j) for j in range(0,12)}
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
    return near+mp.quad(lambda uu: uu**a*mp.log(1-mp.exp(-uu))**b,[delta,1,3,8,20,60,200])
need=[(10-j,j) for j in range(1,11)]
I={ab:I_ab(*ab) for ab in need}
B10=mp.mpf(5)/66
d11 = B10/(10*mp.factorial(10))*sum(mp.binomial(10,j)*(-1)**j*I[(10-j,j)] for j in range(1,11))
print('d11 =', mp.nstr(d11,50))
z=mp.zeta
W11=[z(11),z(2)*z(9),mp.pi**4*z(7),mp.pi**6*z(5),mp.pi**8*z(3),z(2)**2*z(7),z(3)**2*z(5),z(2)**3*z(5),z(2)*z(3)**3]
n11=['z11','z2z9','pi4z7','pi6z5','pi8z3','z2^2z7','z3^2z5','z2^3z5','z2z3^3']
rel=mp.pslq([d11]+W11,tol=mp.mpf('1e-65'),maxcoeff=10**11,maxsteps=10000)
print('d11 PSLQ:', rel)
if rel and rel[0]!=0:
    print(' => d11 = -(', '+'.join(f'{rel[i+1]}*{n11[i]}' for i in range(9)), f')/{rel[0]}')
    val=-sum(mp.mpf(rel[i+1])/rel[0]*b for i,b in enumerate(W11))
    print(' residual:', mp.nstr(abs(d11-val),4))

# (iii) d11 Richardson residual
mp.mp.dps = 105
def A_series(p, tol=mp.mpf('1e-100')):
    p = mp.mpf(p); s = mp.mpf(0); m = 0
    while True:
        m += 1
        t = mp.gamma(1+m/p)**2/(m**2*mp.gamma(1+2*m/p))
        s += t
        if m > 30 and abs(t) < tol: break
    return s
d2 = -mp.zeta(2)/2; d3 = mp.zeta(3)/6
d5 = (mp.zeta(2)*mp.zeta(3)-3*mp.zeta(5))/60
d7 = mp.zeta(7)/14 - mp.zeta(2)*mp.zeta(5)/42 - mp.zeta(2)**2*mp.zeta(3)/140
d9 = (-428400*mp.zeta(9)+136080*mp.zeta(2)*mp.zeta(7)+79*mp.pi**6*mp.zeta(3)+1134*mp.pi**4*mp.zeta(5)-10080*mp.zeta(3)**3)/1814400
xs=[]; ys=[]
for p in list(range(60, 601, 30)):
    e = mp.mpf(1)/p
    D = mp.zeta(2) - A_series(p)
    R = (D - d1*e - d2*e**2 - d3*e**3 - d5*e**5 - d7*e**7 - d9*e**9)/e**11
    xs.append(e); ys.append(R)
deg=3
n=len(xs); M=mp.matrix(n,deg)
for i,x in enumerate(xs):
    for j in range(deg): M[i,j]=x**(2*j)
b=mp.matrix(n,1)
for i,y in enumerate(ys): b[i,0]=y
sol=mp.lu_solve(M.T*M,M.T*b)
print('d11 residual-Richardson:', mp.nstr(sol[0,0],14), ' |diff|', mp.nstr(abs(sol[0,0]-d11),3))
