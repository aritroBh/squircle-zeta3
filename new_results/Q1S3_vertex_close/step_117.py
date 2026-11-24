import pickle, mpmath as mp, json, time
from q1s_exact import femb, BASISEL
t0=time.time()
mp.mp.dps=400
V3=pickle.load(open('q3_accs.pkl','rb'))['V3p']
from collections import Counter
print("word-structure:", Counter((len(w),len(lvm)) for (w,lvm) in V3.W))
print("prod-structure:", Counter(len(lvm) for lvm in V3.P))
# evaluate non-weight-3 parts
import q3_eval, importlib
importlib.reload(q3_eval)
ge=q3_eval.GEval(400, M=96, gln=5, subn=5)
W2=mp.mpc(0)
for (w,lvm),c in V3.W.items():
    if len(w)+len(lvm)==3: continue
    wn=tuple(femb(s,1) for s in w)
    v=ge.word(wn) if len(w)==1 else ge.rest_total(*wn)
    for f in lvm: v*=mp.log(femb(f[1],1))
    W2+=femb(c,1)*v
print("W2 done %.0fs"%(time.time()-t0))
P=mp.mpc(0)
for lvm,c in V3.P.items():
    v=mp.mpc(1)
    for f in lvm: v*=mp.log(femb(f[1],1))
    P+=femb(c,1)*v
print("W2 =", mp.nstr(W2,20))
print("P  =", mp.nstr(P,20))
K=pickle.load(open('K2_V3p_400.pkl','rb'))
Kre=mp.mpf(K['re']); Kim=mp.mpf(K['im'])
K=re(Kre,Kim) if False else mp.mpc(Kre,Kim)
R = K - W2 - P
print("R =", mp.nstr(R,25))
pickle.dump({'re':mp.nstr(mp.re(R),390),'im':mp.nstr(mp.im(R),390)}, open('R_V3p_400.pkl','wb'))
print("time %.0fs"%(time.time()-t0))
