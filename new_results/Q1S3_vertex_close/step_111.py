import pickle, mpmath as mp, json
from q1s_exact import femb, BASISEL
mp.mp.dps=120
sol = pickle.load(open('q3_solve_min.pkl','rb'))['V3p']
sols = pickle.load(open('sunit_li3.pkl','rb'))
X = pickle.load(open('q3_Xli2.pkl','rb'))['X']
def ReLi3(i): return mp.re(mp.polylog(3, femb(sols[i]['z'],1)))
def Li2(i): return mp.polylog(2, femb(sols[i]['z'],1))
def logu(u): return mp.log(femb(BASISEL[u],1))
def Xval(V,W):
    v = logu(V)*logu(W)/2
    for i,c in X[(V,W)]:
        v += mp.mpf(float(c))*Li2(i)
    return v
li3part = sum(mp.mpf(float(c))*2*ReLi3(v[0]) for (k,v),c in sol if k=='li3')
decpart = sum(mp.mpf(float(c))*logu(v[0])*Xval(v[1],v[2]) for (k,v),c in sol if k=='dec')
b3 = "<project>/research/07_push/S3_superellipsoid/claude/"
ref3=mp.mpf(json.load(open(b3+"V3p_1250.json"))['V3p'])
K = ref3 - li3part - decpart
print("li3part:", mp.nstr(li3part,18))
print("decpart:", mp.nstr(decpart,18))
print("K =", mp.nstr(K,25))
