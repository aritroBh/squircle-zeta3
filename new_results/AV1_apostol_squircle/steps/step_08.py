
# step_08.py : d9 PSLQ vs independent weight-9 basis; weight-9 individuals; consolidated table.
import mpmath as mp
mp.mp.dps = 105
import sympy as sp
Kmax = 70
u = sp.symbols('u')
ex = sp.log(sum((-u)**j/sp.factorial(j+1) for j in range(Kmax+3)))
cs = [sp.expand(ex).series(u,0,Kmax+1).removeO().coeff(u,k) for k in range(Kmax+1)]
c = [mp.mpf(x) for x in cs]
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
need9=[(7,1),(6,2),(5,3),(4,4),(3,5),(2,6),(1,7),(0,8)]
I={ab:I_ab(*ab) for ab in need9}
z9,z2z7,p6z3,p4z5,z33 = mp.zeta(9),mp.zeta(2)*mp.zeta(7),mp.pi**6*mp.zeta(3),mp.pi**4*mp.zeta(5),mp.zeta(3)**3
W9=[z9,z2z7,p6z3,p4z5,z33]; names9=['z9','z2z7','pi^6z3','pi^4z5','z3^3']
print('weight-9 individuals:')
for ab in need9:
    print(f' I_{ab}:', mp.pslq([I[ab]]+W9,tol=mp.mpf('1e-40'),maxcoeff=10**8,maxsteps=10000))
d9 = -(-8*I[(7,1)]+28*I[(6,2)]-56*I[(5,3)]+70*I[(4,4)]-56*I[(3,5)]+28*I[(2,6)]-8*I[(1,7)]+I[(0,8)])/9676800
print('d9 =', mp.nstr(d9,40))
r9 = mp.pslq([d9]+W9,tol=mp.mpf('1e-40'),maxcoeff=10**8,maxsteps=10000)
print('d9 PSLQ:', r9, ['1']+names9)
if r9:
    q = r9[0]
    print('i.e. d9 = -(', '+'.join(f'{r9[i+1]}*{names9[i]}' for i in range(5)), f')/{q}')
    val = -sum(sp.Rational(r9[i+1],q)*b for i,b in enumerate(W9))
    print('residual:', mp.nstr(abs(d9-val),5))
# consolidated square-end table
d1 = mp.mpf('1.73300259017088785677509047603076917317936436')
d2 = -mp.zeta(2)/2; d3 = mp.zeta(3)/6
d5 = (mp.zeta(2)*mp.zeta(3)-3*mp.zeta(5))/60
d7 = mp.zeta(7)/14 - mp.zeta(2)*mp.zeta(5)/42 - mp.zeta(2)**2*mp.zeta(3)/140
print()
print('CONSOLIDATED deficit(eps)= zeta2 - A(1/eps):')
print('d1 =', mp.nstr(d1,30), '(OPEN)')
for nm,v in [('d2',d2),('d3',d3),('d5',d5),('d7',d7),('d9',d9)]:
    print(nm,'=',mp.nstr(v,30))
