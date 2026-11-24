import pickle, mpmath as mp, json, time
from q1s_exact import femb, BASISEL
t0=time.time()
DPS=1250
mp.mp.dps=DPS
sols=pickle.load(open('sunit_li3.pkl','rb'))
names=pickle.load(open('q3_names.pkl','rb'))['names']
xnew=pickle.load(open('q3_xnew.pkl','rb'))['xnew']
X=pickle.load(open('q3_Xli2.pkl','rb'))['X']
def ReLi3(i): return mp.re(mp.polylog(3, femb(sols[i]['z'],1)))
def Li2(i): return mp.polylog(2, femb(sols[i]['z'],1))
def logu(u): return mp.log(femb(BASISEL[u],1))
def Xval(V,W):
    v = logu(V)*logu(W)/2
    for i,c in X[(V,W)]: v += mp.mpf(float(c))*Li2(i)
    return v
li3part=mp.mpf(0); decre=mp.mpf(0)
for (k,v),e in zip(names,xnew):
    if e==0: continue
    c=mp.mpf(e)/2
    if k=='li3': li3part += c*2*ReLi3(v[0])
    else: decre += c*mp.re(logu(v[0])*Xval(v[1],v[2]))
print("li3+dec done %.0fs"%(time.time()-t0), flush=True)
V3=pickle.load(open('q3_accs.pkl','rb'))['V3p']
def G2f(r,s):
    if s==0: return mp.polylog(2,1/r)
    if r==0: return -mp.polylog(2,1/s)
    if r==s: return mp.log(1-1/s)**2/2
    ulo=-r/(s-r); uup=(1-r)/(s-r)
    return mp.log((s-r)/s)*(mp.log(uup)-mp.log(ulo)) - (mp.polylog(2,uup)-mp.polylog(2,ulo))
W2=mp.mpc(0)
for (w,lvm),c in V3.W.items():
    if len(w)+len(lvm)==3: continue
    wn=tuple(femb(s,1) for s in w)
    v=mp.log(1-1/wn[0]) if len(w)==1 else G2f(*wn)
    for f in lvm: v*=mp.log(femb(f[1],1))
    W2+=femb(c,1)*v
P=mp.mpc(0)
for lvm,c in V3.P.items():
    v=mp.mpc(1)
    for f in lvm: v*=mp.log(femb(f[1],1))
    P+=femb(c,1)*v
print("W2,P done %.0fs"%(time.time()-t0), flush=True)
b3 = "<project>/research/07_push/S3_superellipsoid/claude/"
ref3=mp.mpf(json.load(open(b3+"V3p_1250.json"))['V3p'])
Rre = ref3 - li3part - decre - mp.re(W2) - mp.re(P)
print("R_real =", mp.nstr(Rre,35))
print("imag check (W2,P):", mp.nstr(mp.im(W2),4), mp.nstr(mp.im(P),4))
pickle.dump({'Rre':mp.nstr(Rre,DPS-20),'li3':mp.nstr(li3part,DPS-20),'dec':mp.nstr(decre,DPS-20),
             'W2':mp.nstr(mp.re(W2),DPS-20),'P':mp.nstr(mp.re(P),DPS-20)}, open('Rreal_V3p_1250.pkl','wb'))
print("time %.0fs"%(time.time()-t0))
