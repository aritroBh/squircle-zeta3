import pickle
from fractions import Fraction as Fr
from q1s_exact import *
V3=pickle.load(open('q3_accs.pkl','rb'))['V3p']
sols=pickle.load(open('sunit_li3.pkl','rb'))
# index sols by (vec_z, vec_1z)
def vkey(v): return tuple(sorted((k,Fr(c)) for k,c in v.items()))
sidx={}
for i,s in enumerate(sols):
    sidx[(vkey(s['vec_z']),vkey(s['vec_1z']))]=i
pairs=set()
for (w,lvm) in V3.W:
    if len(w)==2:
        pairs.add(w)
print("weight-2 word pairs (r,s):", len(pairs))
args=set()
for r,s in pairs:
    if feq(s,ZERO): 
        # G(r,0;1): log-type, Li2 arg = 1/r? G(r,0;1) = int log t/(t-r): involves Li2(1/r)
        args.add(fdiv(ONE,r))
        continue
    u1=fdiv(fsub(ONE,r),fsub(s,r))
    u2=fdiv(fneg(r),fsub(s,r))
    args.add(u1); args.add(u2)
print("cross-ratio args:", len(args))
found=[]
for u in sorted(args,key=repr):
    f=factor_sunit(u)
    if f is None: 
        print("not S-unit:", fshow(u)); continue
    vecu=f[1]
    # also factor 1-u
    f2=factor_sunit(fsub(ONE,u))
    if f2 is None:
        print("1-u not S-unit:", fshow(u)); continue
    i=sidx.get((vkey(vecu),vkey(f2[1])))
    found.append((u,i))
    print("arg", fshow(u), " vec:",vecu, " 1-u:",f2[1], " sol idx:", i)
