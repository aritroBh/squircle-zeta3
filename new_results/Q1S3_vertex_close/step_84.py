import importlib, q3_engine, q1s_exact, q3_stage2
for m in (q1s_exact, q3_engine, q3_stage2): importlib.reload(m)
from q3_stage2 import build_V1, build_V3p
from q1s_exact import *
import pickle
V1=build_V1(); V3=build_V3p()
pickle.dump({'V1':V1,'V3p':V3}, open('q3_accs.pkl','wb'))

# collect all letters appearing in words and all lv arguments
letters=set()
lvs=set()
for acc in (V1,V3):
    for (w,lvm) in acc.W:
        letters.update(w)
        for f in lvm: lvs.add(f[1])
    for lvm in acc.P:
        for f in lvm: lvs.add(f[1])
letters=sorted(letters,key=repr)
print("letters:", [fshow(s) for s in letters])
# check factorizations: differences of consecutive (prev, cur, nxt): cur-prev and nxt-cur... 
# sym_w1(a0,a1,a2) uses (a2-a1)/(a0-a1): need factor_sunit(a2-a1), factor_sunit(a0-a1)
pts = [ZERO]+letters+[ONE]
fail=0
from itertools import permutations
# all pairs (x,y) with x,y in pts x letters... actually needed pairs: (a_{j-1},a_j),(a_{j+1},a_j) for words
pairs=set()
for acc in (V1,V3):
    for (w,lvm) in acc.W:
        p=(ZERO,)+w+(ONE,)
        for j in range(1,len(p)-1):
            pairs.add((p[j-1],p[j])); pairs.add((p[j+1],p[j]))
for x,y in sorted(pairs,key=repr):
    d=fsub(x,y)
    if feq(d,ZERO): continue
    f=factor_sunit(d)
    if f is None:
        fail+=1; print("NOT S-UNIT:", fshow(x)," - ",fshow(y))
    for e in f[1].values():
        pass
print("pairs needing factorization:", len(pairs), " failures:", fail)
# lv arguments
lvfail=0
for x in sorted(lvs,key=repr):
    f=factor_sunit(x)
    if f is None: lvfail+=1; print("lv NOT S-UNIT:", fshow(x))
print("lv args:", len(lvs), " failures:", lvfail)
