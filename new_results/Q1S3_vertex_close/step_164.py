import pickle, mpmath as mp, json, time
from q1s_exact import femb, BASISEL
t0=time.time()
mp.mp.dps=600
sols=pickle.load(open('sunit_li3.pkl','rb'))
names=pickle.load(open('q3_names.pkl','rb'))['names']
x=pickle.load(open('q3_xnew2.pkl','rb'))['xnew2']
X=pickle.load(open('q3_Xli2.pkl','rb'))['X']
def ReLi3(i): return mp.re(mp.polylog(3, femb(sols[i]['z'],1)))
def Li2(i): return mp.polylog(2, femb(sols[i]['z'],1))
def logu(u): return mp.log(femb(BASISEL[u],1))
def Xval(V,W):
    v = logu(V)*logu(W)/2
    for i,c in X[(V,W)]: v += mp.mpf(float(c))*Li2(i)
    return v
decim=mp.mpf(0)
for (k,v),e in zip(names,x):
    if e==0 or k!='dec': continue
    c=mp.mpf(e)/2
    decim += c*mp.im(logu(v[0])*Xval(v[1],v[2]))
print("Im(decpart) =", mp.nstr(decim,25))
r2=mp.sqrt(2); a=1+r2; b=r2-1; L=mp.log(a); lam=mp.log(2); l3=mp.log(3); pi=mp.pi
alpha=mp.atan(1/r2)
for nm,v in [('pi3',pi**3),('pi',pi),('pi*G',pi*mp.catalan),('pi*Cl1',pi*mp.clsin(2,mp.pi/4)),
             ('pi*imLi2p2',pi*mp.im(mp.polylog(2,1-mp.exp(1j*pi/4)))),
             ('pi*reLi2q1',pi*mp.re(mp.polylog(2,r2+1j))),('pi*reLi2p2',pi*mp.re(mp.polylog(2,1-mp.exp(1j*pi/4)))),
             ('pi*lam2',pi*lam**2),('pi*L2',pi*L**2),('pi*lam*L',pi*lam*L),('pi*l3^2',pi*l3**2),
             ('pi*alpha2',pi*alpha**2),('pi*lam*al',pi*lam*alpha),('pi*L*al',pi*L*alpha),('pi*l3*al',pi*l3*alpha),
             ('pi*z2',pi*pi**2/6),('pi*z3',pi*mp.zeta(3))]:
    print("  /%s = %s"%(nm, mp.nstr(decim/v,8)))
print("time %.0fs"%(time.time()-t0))
