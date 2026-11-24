import pickle, mpmath as mp, json, time
from q1s_exact import femb, BASISEL
t0=time.time()
mp.mp.dps=400
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
print("li3part:", mp.nstr(li3part,18), " decre:", mp.nstr(decre,18), " %.0fs"%(time.time()-t0))
b3 = "<project>/research/07_push/S3_superellipsoid/claude/"
ref3=mp.mpf(json.load(open(b3+"V3p_1250.json"))['V3p'])
Kold=pickle.load(open('K2_V3p_400.pkl','rb')); Rold=pickle.load(open('R_V3p_400.pkl','rb'))
W2P = mp.mpf(Kold['re'])-mp.mpf(Rold['re'])
print("W2+P =", mp.nstr(W2P,25))
Rre = ref3 - li3part - decre - W2P
print("R_real =", mp.nstr(Rre,30))
pickle.dump({'Rre':mp.nstr(Rre,390),'li3':mp.nstr(li3part,390),'dec':mp.nstr(decre,390)}, open('Rreal_V3p_400.pkl','wb'))
print("time %.0fs"%(time.time()-t0))
