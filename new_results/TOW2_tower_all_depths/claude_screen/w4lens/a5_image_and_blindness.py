# Items (1)-(3), refined: (A) re-run the exact tower -> kernel linear image for the degree-3 keys at k=4..14
#   from the recorded degree-4 towers (independent of the full reducer), incl. +-i simple pole / const content;
# (B) at k=4, replace the degree-4 towers by their quadratic-in-k extrapolation from k>=5 (all three keys) and
#   check the image kernel equals the exact k=4 kernel  (=> k=4 sits on the quadratic, and is the (2,2) deviation invisible?).
import sys, os, pickle
from fractions import Fraction as Fr
sys.dont_write_bytecode = True
TH="<project>/research/08_theory/WTS_weight_k_structure/theory"
sys.path.insert(0, TH)
from fastred import *
def q4(q): return (q.p0, q.p1, q.p2, q.p3)
def Q(t): return QSI(*t)
def Gm_from_tower(Hd):
    L=[]
    if Hd:
        top=max(Hd); L=[ZERO]*top
        for r,c in Hd.items(): L[r-1]=Q(c)
    H=pf_zero(); H[1]=list(L); H[2]=[e.conj_i() for e in L]
    G=pf_integrate_higher(H); C=pf_eval_at(G,BQ)
    return [G[0],G[1],G[2],-C]
def kernel_from_tower(T,key):
    i,j=key; tot=pf_zero()
    if (i+1,j) in T: tot=pf_add(tot,pf_scale(pf_mul_DL(Gm_from_tower(T[(i+1,j)]),DL1),Fr(-(i+1))))
    if i>j and (i,j+1) in T: tot=pf_add(tot,pf_scale(pf_mul_DL(Gm_from_tower(T[(i,j+1)]),DL2),Fr(-(j+1))))
    if j==i-1 and (i,i) in T: tot=pf_add(tot,pf_scale(pf_fold(pf_mul_DL(Gm_from_tower(T[(i,i)]),DL1)),Fr(-i)))
    return tot
Tw={k:pickle.load(open(os.path.join(TH,f"tower_{k}.pkl"),'rb')) for k in range(3,15)}
Kw={k:pickle.load(open(os.path.join(TH,f"kern_{k}.pkl"),'rb')) for k in range(3,15)}
def realdict(P):
    simple,higher,const=pf_split(P)
    return {(r.p0,r.p1):q4(c) for r,c in simple[0].items() if any(x!=0 for x in q4(c))}, simple, const
print("(A) image of degree-4 towers -> degree-3 kernels")
for k in range(4,15):
    out=[]
    for key in [(3,0),(2,1)]:
        P=kernel_from_tower(Tw[k],key)
        rd,simple,const=realdict(P)
        want={r:c for r,c in Kw[k][key]['real'].items() if any(x!=0 for x in c)}
        has_i=any(not e.is_zero() for e in simple[1]) or any(not e.is_zero() for e in simple[2]) or (not const.is_zero())
        out.append((key, rd==want, 'has_i=%s'%has_i))
    print(" k=%2d"%k, out)
# (B) quadratic extrapolation of degree-4 towers from k>=5 (Lagrange on 5,6,7 per (key,r,comp)); verify on 8..14
def lagr(xs,ys,x):
    s=Fr(0)
    for a,(xi,yi) in enumerate(zip(xs,ys)):
        t=Fr(yi)
        for b,xj in enumerate(xs):
            if b!=a: t*=Fr(x-xj,xi-xj)
        s+=t
    return s
def extrap_tower(k_target, base=(5,6,7), keys=((4,0),(3,1),(2,2))):
    Tx={}
    for key in keys:
        d={}
        for r in range(2,7):
            comps=[]
            for c in range(4):
                ys=[Tw[k][key].get(r,(0,0,0,0))[c] for k in base]
                comps.append(lagr(list(base),ys,k_target))
            d[r]=tuple(comps)
        Tx[key]=d
    return Tx
# check extrapolation reproduces actual towers for k=8..14 (all keys, all r)
badE=0
for k in range(8,15):
    Tx=extrap_tower(k)
    for key in Tx:
        for r in Tx[key]:
            act=tuple(Fr(x) for x in Tw[k][key].get(r,(0,0,0,0)))
            if tuple(Tx[key][r])!=act: badE+=1
print("(B0) quadratic-from-(5,6,7) tower extrapolation mismatches at k=8..14, all 3 degree-4 keys, r=2..6:",badE)
# deviations at k=4
Tx4=extrap_tower(4)
for key in ((4,0),(3,1),(2,2)):
    dev={r:tuple(Fr(Tw[4][key].get(r,(0,0,0,0))[c])-Tx4[key][r][c] for c in range(4)) for r in range(2,7)}
    nzd={r:v for r,v in dev.items() if any(x!=0 for x in v)}
    print(f"(B1) key {key}: deviation (actual k=4 - quadratic extrapolation) nonzero at orders r:",{r:[str(x) for x in v] for r,v in nzd.items()} if nzd else "NONE")
# image of extrapolated towers vs exact k=4 kernels
for key in [(3,0),(2,1)]:
    P=kernel_from_tower(Tx4_ := {kk:{r:tuple(c) for r,c in vv.items()} for kk,vv in Tx4.items()},key)
    # wrap into QSI-friendly (tuples of Fractions are fine for Q())
    rd,simple,const=realdict(P)
    want={r:c for r,c in Kw[4][key]['real'].items() if any(x!=0 for x in c)}
    print(f"(B2) key {key}: image(quadratic-extrapolated degree-4 towers at k=4) == exact k=4 kernel: {rd==want}")
# feed of the (2,2) deviation ALONE into (2,1) kernel at k=4
dev22={r:tuple(Fr(Tw[4][(2,2)].get(r,(0,0,0,0))[c])-Tx4[(2,2)][r][c] for c in range(4)) for r in range(2,7)}
P=kernel_from_tower({(2,2):dev22},(2,1))
rd,simple,const=realdict(P)
print("(B3) feed of the (2,2) k=4 deviation alone into the (2,1) kernel: real residues =",rd,"; +i simple:",[q4(e) for e in simple[1]],"; -i simple:",[q4(e) for e in simple[2]],"; const:",q4(const))
