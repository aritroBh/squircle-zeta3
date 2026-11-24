import pickle, mpmath as mp, itertools, time
from fractions import Fraction as Fr
from q1s_symbol import tadd, tscale, t_conj
from q1s_exact import femb, BASISEL
import sympy as sp
# verify M x = 2t via symbol reconstruction
B4=('a','p2','q1','q2')
names=pickle.load(open('q3_names.pkl','rb'))['names']
sols=pickle.load(open('sunit_li3.pkl','rb'))
xnew=pickle.load(open('q3_xnew.pkl','rb'))['xnew']
S=pickle.load(open('q3_symbols.pkl','rb'))
def li3_sym(pr):
    i,j=pr
    return tadd(sols[i]['S'], t_conj(sols[i]['S']))
def dec_sym(d):
    u,v,w=d
    t={}
    for perm in [(u,v,w),(v,u,w),(v,w,u)]: t[perm]=t.get(perm,Fr(0))+1
    return t
tot={}
for (k,v),e in zip(names,xnew):
    if e==0: continue
    g = li3_sym(v) if k=='li3' else dec_sym(v)
    tot=tadd(tot, tscale(g, Fr(e)))
tgt=S['S3_V3'][0]
ok=all(tot.get(k,Fr(0))==2*tgt.get(k,Fr(0)) for k in set(tot)|set(tgt))
print("M x = 2t exact:", ok)
# evaluate the combination at 120 digits
mp.mp.dps=120
X=pickle.load(open('q3_Xli2.pkl','rb'))['X']
def ReLi3(i): return mp.re(mp.polylog(3, femb(sols[i]['z'],1)))
def Li2(i): return mp.polylog(2, femb(sols[i]['z'],1))
def logu(u): return mp.log(femb(BASISEL[u],1))
def Xval(V,W):
    v = logu(V)*logu(W)/2
    for i,c in X[(V,W)]: v += mp.mpf(float(c))*Li2(i)
    return v
li3part=mp.mpf(0); decpart=mp.mpc(0)
for (k,v),e in zip(names,xnew):
    if e==0: continue
    c=mp.mpf(e)/2
    if k=='li3': li3part += c*2*ReLi3(v[0])
    else: decpart += c*logu(v[0])*Xval(v[1],v[2])
import json
b3 = "<project>/research/07_push/S3_superellipsoid/claude/"
ref3=mp.mpf(json.load(open(b3+"V3p_1250.json"))['V3p'])
K = ref3 - li3part - decpart
print("li3part:", mp.nstr(li3part,15))
print("decpart:", mp.nstr(decpart,15))
print("K =", mp.nstr(K,18))
