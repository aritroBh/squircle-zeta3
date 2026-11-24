import importlib, q3_engine, q1s_exact, q3_eval, q3_stage2
for m in (q1s_exact, q3_engine, q3_eval, q3_stage2): importlib.reload(m)
from q3_stage2 import build_V1, build_V3p
from q3_eval import GEval
from q1s_exact import femb
import mpmath as mp, json, pickle

V1=build_V1(); V3=build_V3p()
pickle.dump({'V1':V1,'V3p':V3}, open('q3_accs.pkl','wb'))

DPS=130
ge=GEval(DPS, M=64, gln=5)
def acc_value(acc, ge):
    tot=mp.mpc(0)
    for (w,lvm),c in acc.W.items():
        wn=tuple(femb(s,1) for s in w)
        v=ge.word(wn)
        for f in lvm: v*=mp.log(femb(f[1],1))
        tot+=femb(c,1)*v
    for lvm,c in acc.P.items():
        v=mp.mpc(1)
        for f in lvm: v*=mp.log(femb(f[1],1))
        tot+=femb(c,1)*v
    tot+=femb(acc.R,1)
    return tot

v1 = acc_value(V1, ge)
v3 = acc_value(V3, ge)
b3 = "<project>/research/07_push/S3_superellipsoid/claude/"
ref1=mp.mpf(json.load(open(b3+"V1_1250.json"))['V1'])
ref3=mp.mpf(json.load(open(b3+"V3p_1250.json"))['V3p'])
print("V1 value  :", mp.nstr(v1,30))
print("V1 imag   :", mp.nstr(mp.im(v1),5))
print("V1 resid  :", mp.nstr(abs(v1-ref1),6), " digits:", mp.nstr(-mp.log10(abs(v1-ref1)),4))
print("V3' value :", mp.nstr(v3,30))
print("V3' imag  :", mp.nstr(mp.im(v3),5))
print("V3' resid :", mp.nstr(abs(v3-ref3),6), " digits:", mp.nstr(-mp.log10(abs(v3-ref3)),4))
