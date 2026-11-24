import pickle, mpmath as mp, time
from q1s_exact import femb, BASISEL
mp.mp.dps=400
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
def re3(i): return mp.re(mp.polylog(3,femb(sols[i]['z'],1)))
k=B[:,1]
v=sum(mp.mpf(int(k[i]))/2*genval[i] for i in range(336) if int(k[i])!=0)
d=v-re3(292)
print("dir1 - reLi3(z292) at 400 digits:", mp.nstr(d,10))
print("|log10|:", mp.nstr(mp.log10(abs(d)) if d!=0 else mp.mpf('-inf'),5))
