import pickle, mpmath as mp, time
from q1s_exact import femb
t0=time.time()
mp.mp.dps=400
sols=pickle.load(open('sunit_li3.pkl','rb'))
r2=mp.sqrt(2); a=1+r2; b=r2-1; L=mp.log(a); lam=mp.log(2); l3=mp.log(3); pi=mp.pi
alpha=mp.atan(1/r2)
z3=mp.zeta(3)
Lchi8=(mp.zeta(3,mp.mpf(1)/8)-mp.zeta(3,mp.mpf(3)/8)-mp.zeta(3,mp.mpf(5)/8)+mp.zeta(3,mp.mpf(7)/8))/512
G=mp.catalan; Cl1=mp.clsin(2,mp.pi/4); z2=pi**2/6
def imli2(z): return mp.im(mp.polylog(2,z))
def reli2(z): return mp.re(mp.polylog(2,z))
def re3(i): return mp.re(mp.polylog(3,femb(sols[i]['z'],1)))
def im3(i): return mp.im(mp.polylog(3,femb(sols[i]['z'],1)))
p2n = 1-mp.exp(1j*pi/4); q1n = r2+1j
base=[('z3',z3),('Lchi8',Lchi8),('pi3',pi**3),('r2*pi3',r2*pi**3),
      ('pi*G',pi*G),('pi*Cl1',pi*Cl1),
      ('z2*lam',z2*lam),('z2*L',z2*L),('z2*l3',z2*l3),('z2*alpha',z2*alpha),
      ('pi*imLi2ib',pi*imli2(1j*b)),('pi*imLi2ia',pi*imli2(1j*a)),
      ('pi*imLi2p2',pi*imli2(p2n)),('pi*imLi2q1',pi*imli2(q1n)),
      ]
for i in [289,292,376,377,399,401,479]:
    base.append(('reLi3(z%d)'%i, re3(i)))
names=[nm for nm,_ in base]; vals=[mp.mpf(v) for _,v in base]
R=pickle.load(open('Rreal_V3p_400.pkl','rb'))
Rre=mp.mpf(R['Rre'])
mp.mp.dps=380
vec=[Rre]+vals
rel=mp.pslq(vec, tol=mp.mpf('1e-340'), maxcoeff=10**20, maxsteps=2000)
ok = rel and rel[0]!=0
print("PSLQ:", rel if ok else "(none)")
if ok:
    res=sum(mp.mpf(r)*v for r,v in zip(rel,vec))
    print("residual:", mp.nstr(abs(res),6))
    print("terms:", [(names[i-1] if i>0 else 'R',r) for i,r in enumerate(rel) if r!=0])
print("time %.0fs"%(time.time()-t0))
