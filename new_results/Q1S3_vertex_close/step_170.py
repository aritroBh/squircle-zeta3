import pickle, mpmath as mp, time
from q1s_exact import femb, BASISEL
t0=time.time()
mp.mp.dps=250
sols=pickle.load(open('sunit_li3.pkl','rb'))
names=pickle.load(open('q3_names.pkl','rb'))['names']
X=pickle.load(open('q3_Xli2.pkl','rb'))['X']
Kred=pickle.load(open('q3_Kred.pkl','rb'))['Kred']
def ReLi3(i): return mp.re(mp.polylog(3, femb(sols[i]['z'],1)))
def Li2(i): return mp.polylog(2, femb(sols[i]['z'],1))
def logu(u): return mp.log(femb(BASISEL[u],1))
def Xval(V,W):
    v = logu(V)*logu(W)/2
    for i,c in X[(V,W)]: v += mp.mpf(float(c))*Li2(i)
    return v
genval=[]
for (k,v) in names:
    if k=='li3': genval.append(2*ReLi3(v[0]))
    else: genval.append(mp.re(logu(v[0])*Xval(v[1],v[2])))
import numpy as np
B=np.array(Kred, dtype=object).T
r2=mp.sqrt(2); a=1+r2; b=r2-1; L=mp.log(a); lam=mp.log(2); l3=mp.log(3); pi=mp.pi
alpha=mp.atan(1/r2)
z3=mp.zeta(3)
Lchi8=(mp.zeta(3,mp.mpf(1)/8)-mp.zeta(3,mp.mpf(3)/8)-mp.zeta(3,mp.mpf(5)/8)+mp.zeta(3,mp.mpf(7)/8))/512
G=mp.catalan; Cl1=mp.clsin(2,mp.pi/4); z2=pi**2/6
def re3(i): return mp.re(mp.polylog(3,femb(sols[i]['z'],1)))
def im3(i): return mp.im(mp.polylog(3,femb(sols[i]['z'],1)))
tors=[243,244,245,246,247,248,249,289,292,376,377,399,401,479]
base=[('z3',z3),('Lchi8',Lchi8),('pi3',pi**3),('r2*pi3',r2*pi**3),('pi*G',pi*G),('pi*Cl1',pi*Cl1),
      ('z2*lam',z2*lam),('z2*L',z2*L),('z2*l3',z2*l3),('z2*alpha',z2*alpha),
      ('Cl3(pi/4)',mp.re(mp.polylog(3,mp.exp(1j*pi/4)))),('Cl3(3pi/4)',mp.re(mp.polylog(3,mp.exp(3j*pi/4))))]
for i in tors:
    base.append(('reLi3(z%d)'%i, re3(i)))
    base.append(('imLi3(z%d)'%i, im3(i)))
base=[(nm,mp.mpf(v)) for nm,v in base if abs(mp.mpf(v))>mp.mpf('1e-230')]
names2=[nm for nm,_ in base]; vals=[v for _,v in base]
for j in [0,1,2,3,4,5]:
    k=B[:,j]
    v=sum(mp.mpf(int(k[i]))/2*genval[i] for i in range(336) if int(k[i])!=0)
    vec=[v]+vals
    rel=mp.pslq(vec, tol=mp.mpf('1e-220'), maxcoeff=10**10, maxsteps=2000)
    ok = rel and rel[0]!=0
    print("dir",j, "value", mp.nstr(v,14), ":", "FOUND" if ok else "(none)")
    if ok:
        res=sum(mp.mpf(r)*vv for r,vv in zip(rel,vec))
        print("   residual:", mp.nstr(abs(res),6), " terms:", [(names2[i-1] if i>0 else 'K',r) for i,r in enumerate(rel) if r!=0])
print("time %.0fs"%(time.time()-t0))
