# Items (2),(3): exact residue-by-residue polynomial-in-k analysis of the degree-3 kernels, all levels 3..21.
import pickle, os, sys
from fractions import Fraction as Fr
sys.dont_write_bytecode = True
TH="<project>/research/08_theory/WTS_weight_k_structure/theory"
RD="<project>/research/08_theory/TWOEND1_weight_poles/red"
def path(k):
    for D in (TH,RD):
        p=os.path.join(D,f"kern_{k}.pkl")
        if os.path.exists(p): return p
K={k:pickle.load(open(path(k),'rb')) for k in range(3,22)}
def lagr(xs, ys, x):
    # exact Lagrange interpolation
    s=Fr(0)
    for i,(xi,yi) in enumerate(zip(xs,ys)):
        t=Fr(yi)
        for j,xj in enumerate(xs):
            if j!=i: t*=Fr(x-xj, xi-xj)
        s+=t
    return s
def third_diff_zero(vals):
    # vals list at consecutive k; return True iff 3rd differences vanish
    d=vals
    for _ in range(3): d=[d[i+1]-d[i] for i in range(len(d)-1)]
    return all(x==0 for x in d)
def second_diff_nonzero(vals):
    d=vals
    for _ in range(2): d=[d[i+1]-d[i] for i in range(len(d)-1)]
    return any(x!=0 for x in d)
for key in [(3,0),(2,1)]:
    print("===== key",key)
    poles=set()
    for k in range(3,22): poles|=set(K[k][key]['real'].keys())
    print("poles (p0,p1)=p0+p1*sqrt2 appearing at any k in 3..21:",sorted(poles))
    for pole in sorted(poles):
        for comp in range(4):
            val={k: Fr(K[k][key]['real'].get(pole,(0,0,0,0))[comp]) for k in range(3,22)}
            if all(v==0 for v in val.values()): continue
            # quadratic through k=5,6,7
            xs=[5,6,7]; ys=[val[k] for k in xs]
            bad_hi=[k for k in range(8,22) if lagr(xs,ys,k)!=val[k]]
            e4=val[4]-lagr(xs,ys,4); e3=val[3]-lagr(xs,ys,3)
            # also: 3rd differences on k=5..21 and on k=4..21
            v5=[val[k] for k in range(5,22)]; v4=[val[k] for k in range(4,22)]
            print(f" pole={pole} comp={comp}: quad(5,6,7) matches k=8..21: {not bad_hi}; "
                  f"k=4 residual (actual - quad)={e4}; k=3 residual={e3}; 3rd-diff-zero on 5..21: {third_diff_zero(v5)}; on 4..21: {third_diff_zero(v4)}; true degree 2: {second_diff_nonzero(v5)}")
    print(" values k=4..7 for reference (pole,comp)->[k4,k5,k6,k7]:")
    for pole in sorted(poles):
        for comp in range(4):
            vv=[Fr(K[k][key]['real'].get(pole,(0,0,0,0))[comp]) for k in (3,4,5,6,7)]
            if any(vv): print("   ",pole,comp,[str(x) for x in vv])
# any pole at +b (-1,1)?  at all degrees / levels
cnt=0
for k in range(1,22):
    for key,e in K.get(k,{}).items():
        if (Fr(-1),Fr(1)) in e['real']:
            if any(x!=0 for x in e['real'][(Fr(-1),Fr(1))]): cnt+=1; print("POLE AT +b:",k,key)
print("kernels (k=3..21, all keys) with nonzero residue at +b:",cnt)
# reality: components 2,3 zero for all residues of degree-3 kernels
bad=0
for k in range(3,22):
    for key in [(3,0),(2,1)]:
        for pole,c in K[k][key]['real'].items():
            if c[2]!=0 or c[3]!=0: bad+=1
print("non-real residue components in degree-3 kernels:",bad)
