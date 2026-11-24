import importlib, q3_engine, q1s_exact, q3_eval, q3_stage2
for m in (q1s_exact, q3_engine, q3_eval, q3_stage2): importlib.reload(m)
from q3_stage2 import build_V1, build_V3p
from q3_eval import GEval
from q1s_exact import femb
import mpmath as mp, json, pickle, time

V1=build_V1(); V3=build_V3p()
pickle.dump({'V1':V1,'V3p':V3}, open('q3_accs.pkl','wb'))
DPS=130
t0=time.time()
ge=GEval(DPS, M=64, gln=5, subn=5)
allw=set()
for acc in (V1,V3):
    for (w,lvm) in acc.W:
        if len(w)>=2:
            allw.add(tuple(femb(s,1) for s in w))
w3 = ge.word_batch([w for w in allw if len(w)==3])
w2 = {w: ge.rest_total(*w) for w in allw if len(w)==2}
print("word eval time:", round(time.time()-t0,1),"s;  unique words:", len(allw))
def acc_value(acc):
    tot=mp.mpc(0)
    for (w,lvm),c in acc.W.items():
        wn=tuple(femb(s,1) for s in w)
        v=ge.word(wn) if len(w)==1 else (w2[wn] if len(w)==2 else w3[wn])
        for f in lvm: v*=mp.log(femb(f[1],1))
        tot+=femb(c,1)*v
    for lvm,c in acc.P.items():
        v=mp.mpc(1)
        for f in lvm: v*=mp.log(femb(f[1],1))
        tot+=femb(c,1)*v
    tot+=femb(acc.R,1)
    return tot
v1=acc_value(V1); v3=acc_value(V3)
b3 = "<project>/research/07_push/S3_superellipsoid/claude/"
ref1=mp.mpf(json.load(open(b3+"V1_1250.json"))['V1'])
ref3=mp.mpf(json.load(open(b3+"V3p_1250.json"))['V3p'])
print("V1  =", mp.nstr(v1,25))
print("V1  resid:", mp.nstr(abs(v1-ref1),6), " digits:", mp.nstr(-mp.log10(abs(v1-ref1)),5))
print("V3' =", mp.nstr(v3,25))
print("V3' resid:", mp.nstr(abs(v3-ref3),6), " digits:", mp.nstr(-mp.log10(abs(v3-ref3)),5))
