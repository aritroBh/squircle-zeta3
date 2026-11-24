import pickle, mpmath as mp, json, time
from fractions import Fraction as Fr
from q1s_symbol import tadd, tscale, t_conj
from q1s_exact import femb, BASISEL
# verify symbol identity for xnew2
B4=('a','p2','q1','q2')
names=pickle.load(open('q3_names.pkl','rb'))['names']
sols=pickle.load(open('sunit_li3.pkl','rb'))
x=pickle.load(open('q3_xnew2.pkl','rb'))['xnew2']
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
for (k,v),e in zip(names,x):
    if e==0: continue
    g = li3_sym(v) if k=='li3' else dec_sym(v)
    tot=tadd(tot, tscale(g, Fr(e)))
tgt=S['S3_V3'][0]
ok=all(tot.get(k,Fr(0))==2*tgt.get(k,Fr(0)) for k in set(tot)|set(tgt))
print("M x = 2t exact:", ok)
# compute R at 1250 digits
t0=time.time()
mp.mp.dps=1250
X=pickle.load(open('q3_Xli2.pkl','rb'))['X']
def ReLi3(i): return mp.re(mp.polylog(3, femb(sols[i]['z'],1)))
def Li2(i): return mp.polylog(2, femb(sols[i]['z'],1))
def logu(u): return mp.log(femb(BASISEL[u],1))
def Xval(V,W):
    v = logu(V)*logu(W)/2
    for i,c in X[(V,W)]: v += mp.mpf(float(c))*Li2(i)
    return v
li3part=mp.mpf(0); decre=mp.mpf(0)
for (k,v),e in zip(names,x):
    if e==0: continue
    c=mp.mpf(e)/2
    if k=='li3': li3part += c*2*ReLi3(v[0])
    else: decre += c*mp.re(logu(v[0])*Xval(v[1],v[2]))
print("parts %.0fs"%(time.time()-t0), flush=True)
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
ref3=mp.mpf(json.load(open("<project>/research/07_push/S3_superellipsoid/claude/V3p_1250.json"))['V3p'])
Rre = ref3 - li3part - decre - mp.re(W2) - mp.re(P)
print("R_real =", mp.nstr(Rre,40))
pickle.dump({'Rre':mp.nstr(Rre,1230)}, open('Rreal2_V3p_1250.pkl','wb'))
print("time %.0fs"%(time.time()-t0))
