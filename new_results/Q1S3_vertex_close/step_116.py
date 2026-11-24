import pickle, mpmath as mp
from q1s_exact import femb
mp.mp.dps=120
K=pickle.load(open('K2_V3p_400.pkl','rb'))
Kre=mp.mpf(K['re']); Kim=mp.mpf(K['im'])
sols=pickle.load(open('sunit_li3.pkl','rb'))
r2=mp.sqrt(2); a=1+r2; b=r2-1; L=mp.log(a); lam=mp.log(2); l3=mp.log(3); pi=mp.pi
alpha=mp.atan(1/r2)
z3=mp.zeta(3)
Lchi8=(mp.zeta(3,mp.mpf(1)/8)-mp.zeta(3,mp.mpf(3)/8)-mp.zeta(3,mp.mpf(5)/8)+mp.zeta(3,mp.mpf(7)/8))/512
G=mp.catalan; Cl1=mp.clsin(2,mp.pi/4); z2=pi**2/6
p2n = 1-mp.exp(1j*pi/4); q1n = r2+1j
Ti2 = lambda x: (mp.polylog(2,1j*x)-mp.polylog(2,-1j*x))/(2j)
cands={'z3':z3,'Lchi8':Lchi8,'pi3':pi**3,'pi*G':pi*G,'pi*Cl1':pi*Cl1,
 'pi*reLi2q1':pi*mp.re(mp.polylog(2,q1n)),'pi*imLi2q1':pi*mp.im(mp.polylog(2,q1n)),
 'pi*Ti2(r2)':pi*Ti2(r2),'pi*Ti2(1/r2)':pi*Ti2(1/r2),
 'pi*Cl2(2a)':pi*mp.clsin(2,2*alpha),'pi*Cl2(pi-2a)':pi*mp.clsin(2,pi-2*alpha),
 'pi*L2':pi*L**2,'pi*lam2':pi*lam**2,'pi*lam*L':pi*lam*L,'pi*alpha2':pi*alpha**2,
 'pi*imLi2ib':pi*mp.im(mp.polylog(2,1j*b)),'pi*imLi2ia':pi*mp.im(mp.polylog(2,1j*a)),
 'pi*imLi2p2':pi*mp.im(mp.polylog(2,p2n)),'z2*alpha':z2*alpha,'pi*z3':pi*z3,
 'Cl3(pi/4)':mp.re(mp.polylog(3,mp.exp(1j*pi/4))),'pi*reLi2p2':pi*mp.re(mp.polylog(2,p2n))}
print("Im K =", mp.nstr(Kim,15))
for nm,v in cands.items():
    print("  Im/%s = %s"%(nm, mp.nstr(Kim/v,8)))
