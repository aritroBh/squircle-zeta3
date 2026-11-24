import pickle, mpmath as mp, time
from q1s_exact import femb
t0=time.time()
mp.mp.dps=1150
sols=pickle.load(open('sunit_li3.pkl','rb'))
r2=mp.sqrt(2); a=1+r2; b=r2-1; L=mp.log(a); lam=mp.log(2); l3=mp.log(3); pi=mp.pi
alpha=mp.atan(1/r2)
z3=mp.zeta(3)
Lchi8=(mp.zeta(3,mp.mpf(1)/8)-mp.zeta(3,mp.mpf(3)/8)-mp.zeta(3,mp.mpf(5)/8)+mp.zeta(3,mp.mpf(7)/8))/512
G=mp.catalan; Cl1=mp.clsin(2,mp.pi/4); z2=pi**2/6
def re3(i): return mp.re(mp.polylog(3,femb(sols[i]['z'],1)))
def im3(i): return mp.im(mp.polylog(3,femb(sols[i]['z'],1)))
R=pickle.load(open('Rreal2_V3p_1250.pkl','rb'))
Rre=mp.mpf(R['Rre'])
tors=[243,244,245,246,247,248,249,289,292,376,377,399,401,479]
base=[('z3',z3),('Lchi8',Lchi8),('pi3',pi**3),('r2*pi3',r2*pi**3),
      ('Cl3(pi/4)',mp.re(mp.polylog(3,mp.exp(1j*pi/4)))),('Cl3(3pi/4)',mp.re(mp.polylog(3,mp.exp(3j*pi/4))))]
for i in tors:
    base.append(('reLi3(z%d)'%i, re3(i)))
    base.append(('imLi3(z%d)'%i, im3(i)))
base=[(nm,mp.mpf(v)) for nm,v in base if abs(mp.mpf(v))>mp.mpf('1e-1050')]
names2=[nm for nm,_ in base]; vals=[v for _,v in base]
print("basis:", len(vals), "-> max precision needed ~", 10*len(vals))
vec=[Rre]+vals
rel=mp.pslq(vec, tol=mp.mpf('1e-1000'), maxcoeff=10**10, maxsteps=3000)
ok = rel and rel[0]!=0
print("PSLQ:", "FOUND" if ok else "(none)")
if ok:
    res=sum(mp.mpf(r)*v for r,v in zip(rel,vec))
    print("residual:", mp.nstr(abs(res),6))
    print("terms:", [(names2[i-1] if i>0 else 'R',r) for i,r in enumerate(rel) if r!=0])
print("time %.0fs"%(time.time()-t0))
