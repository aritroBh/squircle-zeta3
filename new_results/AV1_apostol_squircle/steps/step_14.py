
# step_14.py : zeta(3) content of c_3; small-pool PSLQ for c_1, c_2-residual, W=c_3+2*zeta3*Sm.
import mpmath as mp
mp.mp.dps = 100
def f(m): return mp.gamma(m+1)**2/mp.gamma(2*m+1)
N=900
H1=[mp.mpf(0)];H2=[mp.mpf(0)];H3=[mp.mpf(0)]
for m in range(1,2*N+1):
    H1.append(H1[-1]+1/mp.mpf(m));H2.append(H2[-1]+1/mp.mpf(m)**2);H3.append(H3[-1]+1/mp.mpf(m)**3)
def gder(m,j):
    if j==1: return 2*(H1[m]-mp.euler) - 2*(H1[2*m]-mp.euler)
    if j==2: return -2*mp.zeta(2)-2*H2[m]+4*H2[2*m]
    if j==3: return 12*mp.zeta(3)+4*H3[m]-16*H3[2*m]
c1=c2=c3=mp.mpf(0); Sm=S0=mp.mpf(0)
for m in range(1,N+1):
    fm=f(m)
    if fm<mp.mpf('1e-115') and m>60: break
    S0+=fm; Sm+=m*fm
    g1,g2,g3 = gder(m,1),gder(m,2),gder(m,3)
    c1 += -1*(fm/m)*g1
    c2 += mp.mpf('0.5')*fm*(g2+g1**2)
    c3 += -(mp.mpf(m)*fm/6)*(g3+3*g1*g2+g1**3)
rt3=mp.sqrt(3); Cl=mp.clsin(2,mp.pi/3)
print('c1 =', mp.nstr(c1,30))
print('c2 =', mp.nstr(c2,30))
print('c3 =', mp.nstr(c3,30))
print('Sm = (6+2pi/rt3)/9 check:', mp.nstr(abs(Sm-(6+2*mp.pi/rt3)/9),3))
print('S0 = (3+2pi/rt3)/9 check:', mp.nstr(abs(S0-(3+2*mp.pi/rt3)/9),3))
W = c3 + 2*mp.zeta(3)*Sm
print('W = c3 + 2*zeta3*Sm =', mp.nstr(W,30))
pool=[mp.mpf(1),mp.pi/rt3,mp.log(2),mp.log(3),mp.zeta(2),mp.zeta(3),Cl,Cl*mp.pi,mp.zeta(3)*mp.pi]
bn=['1','pi/rt3','log2','log3','z2','z3','Cl','pi*Cl','z3*pi']
for nm,v in [('W',W),('c1',c1),('c3',c3)]:
    print(nm,'PSLQ:', mp.pslq([v]+pool,tol=mp.mpf('1e-38'),maxcoeff=10**6,maxsteps=10000))
# small pools for c1
for pool2,nm2 in [([mp.mpf(1),mp.pi/rt3,mp.log(2)],'1,pi/rt3,log2'),
                  ([mp.mpf(1),mp.pi/rt3],'*1,pi/rt3'),
                  ([mp.log(2),mp.pi/rt3],'*log2,pi/rt3'),
                  ([mp.mpf(1),mp.log(2)],'*1,log2')]:
    print('c1 PSLQ vs',nm2,':', mp.pslq([c1]+pool2,tol=mp.mpf('1e-42'),maxcoeff=10**9,maxsteps=10000))
# c2 residual: c2 + zeta2*S0 (the provable zeta2 part), then PSLQ
R2 = c2 + mp.zeta(2)*S0
print('R2 = c2+zeta2*S0 =', mp.nstr(R2,30))
print('R2 PSLQ:', mp.pslq([R2]+pool,tol=mp.mpf('1e-38'),maxcoeff=10**6,maxsteps=10000))
print('R2 PSLQ small {z2,pi/rt3,1}:', mp.pslq([R2,mp.zeta(2),mp.pi/rt3,mp.mpf(1)],tol=mp.mpf('1e-42'),maxcoeff=10**9,maxsteps=10000))
