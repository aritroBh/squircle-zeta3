import pickle, mpmath as mp, time
from q1s_exact import femb, BASISEL
t0=time.time()
mp.mp.dps=200
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
vals=[]
for j in [0,1,2,3,4,5,10,50,100,200,280]:
    k=B[:,j]
    v=sum(mp.mpf(int(k[i]))/2*genval[i] for i in range(336) if int(k[i])!=0)
    vals.append(v)
    print("kernel dir %3d: value = %s"%(j, mp.nstr(v,20)))
print("time %.0fs"%(time.time()-t0))
