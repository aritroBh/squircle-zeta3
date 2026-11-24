import pickle
from fractions import Fraction as Fr
from q1s_exact import *
def fdiv(x,y): return fmul(x,finv(y))
V3=pickle.load(open('q3_accs.pkl','rb'))['V3p']
sols=pickle.load(open('sunit_li3.pkl','rb'))
def vkey(v): return tuple(sorted((k,Fr(c)) for k,c in v.items()))
sidx={}
for i,s in enumerate(sols):
    sidx[(vkey(s['vec_z']),vkey(s['vec_1z']))]=i
pairs=set()
for (w,lvm) in V3.W:
    if len(w)==2:
        pairs.add(w)
args=set()
for r,s in pairs:
    if feq(s,ZERO): 
        args.add(fdiv(ONE,r)); continue
    args.add(fdiv(fsub(ONE,r),fsub(s,r)))
    args.add(fdiv(fneg(r),fsub(s,r)))
print("pairs:",len(pairs)," args:",len(args))
found=[]
for u in sorted(args,key=repr):
    f=factor_sunit(u)
    f2=factor_sunit(fsub(ONE,u)) if f else None
    if f is None or f2 is None:
        print("NOT S-UNIT:", fshow(u), f, f2); continue
    i=sidx.get((vkey(f[1]),vkey(f2[1])))
    found.append(i)
print("args not in 493:", [fshow(u) for u,i in zip(sorted(args,key=repr),found) if i is None] if any(i is None for _,i in found) else "none")
print("indices found:", sorted(set(i for i in found if i is not None)))
