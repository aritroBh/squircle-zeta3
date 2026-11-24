import pickle, mpmath as mp, time
from q1s_exact import femb, BASISEL
t0=time.time()
mp.mp.dps=120
sols=pickle.load(open('sunit_li3.pkl','rb'))
names=pickle.load(open('q3_names.pkl','rb'))['names']
X=pickle.load(open('q3_Xli2.pkl','rb'))['X']
Kred=pickle.load(open('q3_Kred.pkl','rb'))['Kred']
xnew=pickle.load(open('q3_xnew.pkl','rb'))['xnew']
def ReLi3(i): return mp.re(mp.polylog(3, femb(sols[i]['z'],1)))
def Li2(i): return mp.polylog(2, femb(sols[i]['z'],1))
def logu(u): return mp.log(femb(BASISEL[u],1))
def Xval(V,W):
    v = logu(V)*logu(W)/2
    for i,c in X[(V,W)]: v += mp.mpf(float(c))*Li2(i)
    return v
# value vector of each generator (li3: 2ReLi3, dec: Re(log u X))
genval=[]
for (k,v) in names:
    if k=='li3': genval.append(2*ReLi3(v[0]))
    else: genval.append(mp.re(logu(v[0])*Xval(v[1],v[2])))
import json
b3 = "<project>/research/07_push/S3_superellipsoid/claude/"
ref3=mp.mpf(json.load(open(b3+"V3p_1250.json"))['V3p'])
W2P = mp.mpf('1.495981428555097652986385')
import numpy as np
B=np.array(Kred, dtype=object).T   # 336x281 kernel
x=np.array(xnew, dtype=object)
def Rval(x):
    s=ref3
    for i in range(336):
        e=int(x[i])
        if e: s -= mp.mpf(e)/2*genval[i]
    return s - W2P
cur=Rval(x); curn=abs(cur)
print("start |R|:", mp.nstr(curn,8), " %.0fs"%(time.time()-t0))
# greedy: move along kernel vectors to reduce |R|
improved=True; it=0
while improved and it<300:
    improved=False; it+=1
    for j in range(281):
        for s in (1,-1):
            y=x+s*B[:,j]
            v=Rval(y)
            if abs(v) < curn*mp.mpf('0.999'):
                x=y; cur=v; curn=abs(v); improved=True
print("after greedy: |R| =", mp.nstr(curn,10), " R =", mp.nstr(cur,18))
xl=[int(e) for e in x]
print("max|x|:", max(abs(e) for e in xl), " norm2:", sum(e*e for e in xl))
pickle.dump({'xnew2':xl}, open('q3_xnew2.pkl','wb'))
print("time %.0fs"%(time.time()-t0))
