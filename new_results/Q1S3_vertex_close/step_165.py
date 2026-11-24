import pickle, mpmath as mp, time
from q1s_exact import femb
t0=time.time()
mp.mp.dps=600
sols=pickle.load(open('sunit_li3.pkl','rb'))
names=pickle.load(open('q3_names.pkl','rb'))['names']
x=pickle.load(open('q3_xnew2.pkl','rb'))['xnew2']
X=pickle.load(open('q3_Xli2.pkl','rb'))['X']
def Li2(i): return mp.polylog(2, femb(sols[i]['z'],1))
def logu(u): return mp.log(femb(BASISEL[u],1))
def Xval(V,W):
    v = logu(V)*logu(W)/2
    for i,c in X[(V,W)]: v += mp.mpf(float(c))*Li2(i)
    return v
decim=mp.mpf(0)
for (k,v),e in zip(names,x):
    if e==0 or k!='dec': continue
    decim += mp.mpf(e)/2*mp.im(logu(v[0])*Xval(v[1],v[2]))
r2=mp.sqrt(2); a=1+r2; b=r2-1; L=mp.log(a); lam=mp.log(2); l3=mp.log(3); pi=mp.pi
alpha=mp.atan(1/r2)
z3=mp.zeta(3)
Lchi8=(mp.zeta(3,mp.mpf(1)/8)-mp.zeta(3,mp.mpf(3)/8)-mp.zeta(3,mp.mpf(5)/8)+mp.zeta(3,mp.mpf(7)/8))/512
G=mp.catalan; Cl1=mp.clsin(2,mp.pi/4); z2=pi**2/6
def imli2(z): return mp.im(mp.polylog(2,z))
def reli2(z): return mp.re(mp.polylog(2,z))
p2n = 1-mp.exp(1j*pi/4); q1n = r2+1j
base=[('pi3',pi**3),('z3',z3),('Lchi8',Lchi8),('pi*G',pi*G),('pi*Cl1',pi*Cl1),
      ('z2*lam',z2*lam),('z2*L',z2*L),('z2*l3',z2*l3),('z2*alpha',z2*alpha),
      ('pi*lam2',pi*lam**2),('pi*L2',pi*L**2),('pi*lam*L',pi*lam*L),('pi*l3^2',pi*l3**2),
      ('pi*alpha2',pi*alpha**2),('pi*lam*al',pi*lam*alpha),('pi*L*al',pi*L*alpha),('pi*l3*al',pi*l3*alpha),
      ('pi*reLi2p2',pi*reli2(p2n)),('pi*reLi2q1',pi*reli2(q1n)),
      ('pi*reLi2a',pi*reli2(a)),('pi*reLi2ib',pi*reli2(1j*b)),('pi*reLi2ia',pi*reli2(1j*a)),
      ]
zused = sorted(set(v[0] for (k,v),e in zip(names,x) if k=='li3' and e!=0)+[i for pr in X.values() for i,c in pr])
for i in zused:
    base.append(('pi*reLi2(z%d)'%i, pi*reli2(femb(sols[i]['z'],1))))
base=[(nm,mp.mpf(v)) for nm,v in base if abs(mp.mpf(v))>mp.mpf('1e-550')]
names2=[nm for nm,_ in base]; vals=[v for _,v in base]
print("basis:", len(vals))
vec=[decim]+vals
rel=mp.pslq(vec, tol=mp.mpf('1e-520'), maxcoeff=10**30, maxsteps=3000)
ok = rel and rel[0]!=0
print("PSLQ:", "FOUND" if ok else "(none)")
if ok:
    res=sum(mp.mpf(r)*v for r,v in zip(rel,vec))
    print("residual:", mp.nstr(abs(res),6))
    print("terms:", [(names2[i-1] if i>0 else 'Im',r) for i,r in enumerate(rel) if r!=0])
print("time %.0fs"%(time.time()-t0))
