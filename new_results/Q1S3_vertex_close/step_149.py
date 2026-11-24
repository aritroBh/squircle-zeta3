import pickle, mpmath as mp, time
from q1s_exact import femb, fsub, fmul, finv, ONE
t0=time.time()
mp.mp.dps=400
sols=pickle.load(open('sunit_li3.pkl','rb'))
names=pickle.load(open('q3_names.pkl','rb'))['names']
xnew=pickle.load(open('q3_xnew.pkl','rb'))['xnew']
X=pickle.load(open('q3_Xli2.pkl','rb'))['X']
r2=mp.sqrt(2); a=1+r2; b=r2-1; L=mp.log(a); lam=mp.log(2); l3=mp.log(3); pi=mp.pi
alpha=mp.atan(1/r2)
z3=mp.zeta(3)
Lchi8=(mp.zeta(3,mp.mpf(1)/8)-mp.zeta(3,mp.mpf(3)/8)-mp.zeta(3,mp.mpf(5)/8)+mp.zeta(3,mp.mpf(7)/8))/512
G=mp.catalan; Cl1=mp.clsin(2,mp.pi/4); z2=pi**2/6
def imli2(z): return mp.im(mp.polylog(2,z))
def reli2(z): return mp.re(mp.polylog(2,z))
def Ti2(x): return mp.re((mp.polylog(2,1j*x)-mp.polylog(2,-1j*x))/(2j))
def re3(i): return mp.re(mp.polylog(3,femb(sols[i]['z'],1)))
def im3(i): return mp.im(mp.polylog(3,femb(sols[i]['z'],1)))
p2n = 1-mp.exp(1j*pi/4); q1n = r2+1j
zused = sorted(set([v[0] for (k,v),e in zip(names,xnew) if k=='li3' and e!=0]+[i for pr in X.values() for i,c in pr]))
cr = [13, 37, 38, 93, 144, 150, 151, 184, 185, 223, 227, 238, 242, 282, 286, 291, 345, 347, 348, 370, 372, 380]
base=[('z3',z3),('Lchi8',Lchi8),('pi3',pi**3),('r2*pi3',r2*pi**3),
      ('Cl3(pi/4)',mp.re(mp.polylog(3,mp.exp(1j*pi/4)))),('Cl3(3pi/4)',mp.re(mp.polylog(3,mp.exp(3j*pi/4)))),
      ('pi*G',pi*G),('pi*Cl1',pi*Cl1),
      ('z2*lam',z2*lam),('z2*L',z2*L),('z2*l3',z2*l3),('z2*alpha',z2*alpha),
      ('pi*lam2',pi*lam**2),('pi*L2',pi*L**2),('pi*lam*L',pi*lam*L),('pi*alpha2',pi*alpha**2),
      ('pi*lam*al',pi*lam*alpha),('pi*L*al',pi*L*alpha),('pi*l3*al',pi*l3*alpha),
      ('pi*imLi2ib',pi*imli2(1j*b)),('pi*imLi2ia',pi*imli2(1j*a)),
      ('pi*imLi2p2',pi*imli2(p2n)),('pi*imLi2q1',pi*imli2(q1n)),
      ('pi*reLi2p2',pi*reli2(p2n)),('pi*reLi2q1',pi*reli2(q1n)),
      ('pi*reLi2a',pi*reli2(a)),('pi*reLi2ib',pi*reli2(1j*b)),('pi*reLi2ia',pi*reli2(1j*a)),
      ('pi*Ti2(r2)',pi*Ti2(r2)),('pi*Cl2(2a)',pi*mp.clsin(2,2*alpha)),('pi*Cl2(pi-2a)',pi*mp.clsin(2,pi-2*alpha)),
      ]
tors=[243,244,245,246,247,248,249,289,292,376,377,399,401,479]
for i in tors:
    base.append(('reLi3(z%d)'%i, re3(i)))
    base.append(('imLi3(z%d)'%i, im3(i)))
# the 2 unsolved cross-ratio args as raw Li2 (compute their embeddings from F-tuples)
from fractions import Fraction as Fr
u1=(Fr(1,2),Fr(-1,4),Fr(0),Fr(1,4)); u2=(Fr(1,2),Fr(1,4),Fr(0),Fr(-1,4))
base.append(('pi*reLi2(u1)', pi*reli2(femb(u1,1))))
base.append(('pi*imLi2(u1)', pi*imli2(femb(u1,1))))
base.append(('pi*reLi2(u2)', pi*reli2(femb(u2,1))))
base.append(('pi*imLi2(u2)', pi*imli2(femb(u2,1))))
for i in sorted(set(zused+cr)):
    base.append(('pi*reLi2(z%d)'%i, pi*reli2(femb(sols[i]['z'],1))))
    base.append(('pi*imLi2(z%d)'%i, pi*imli2(femb(sols[i]['z'],1))))
base=[(nm,mp.mpf(v)) for nm,v in base if abs(mp.mpf(v))>mp.mpf('1e-300')]
names2=[nm for nm,_ in base]; vals=[v for _,v in base]
print("basis:", len(vals), " %.0fs"%(time.time()-t0))
R=pickle.load(open('Rreal_V3p_400.pkl','rb'))
Rre=mp.mpf(R['Rre'])
mp.mp.dps=380
vec=[Rre]+vals
rel=mp.pslq(vec, tol=mp.mpf('1e-290'), maxcoeff=10**22, maxsteps=3000)
ok = rel and rel[0]!=0
print("PSLQ:", "FOUND" if ok else "(none)")
if ok:
    res=sum(mp.mpf(r)*v for r,v in zip(rel,vec))
    print("residual:", mp.nstr(abs(res),6))
    print("terms:", [(names2[i-1] if i>0 else 'R',r) for i,r in enumerate(rel) if r!=0])
print("time %.0fs"%(time.time()-t0))
