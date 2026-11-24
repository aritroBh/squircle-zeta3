# Item (1): degree-4 tower coefficients vs k. Exact Fractions. Check polynomial degree floor(delta/2), and from which k.
import pickle, os, sys
from fractions import Fraction as Fr
sys.dont_write_bytecode = True
TH="<project>/research/08_theory/WTS_weight_k_structure/theory"
T={}
for k in range(2,15):
    p=os.path.join(TH,f"tower_{k}.pkl")
    if os.path.exists(p): T[k]=pickle.load(open(p,'rb'))
print("tower levels:",sorted(T))
print("keys at k=6:",sorted(T[6].keys()))
def diffs(vals):
    out=[vals]
    while len(out[-1])>1:
        out.append([out[-1][i+1]-out[-1][i] for i in range(len(out[-1])-1)])
    return out
def polydeg(ks, vals):
    # smallest m s.t. m-th forward difference identically 0 (needs >= m+2 points to claim)
    D=diffs(vals)
    for m,row in enumerate(D):
        if all(x==0 for x in row): return m-1 if m>0 else -1 , len(vals)
    return None, len(vals)
n=4
for key in [(4,0),(3,1),(2,2)]:
    print("=== key",key)
    for r in range(2,2*n-1):
        delta=2*n-2-r
        for comp in range(4):
            ks=[k for k in sorted(T) if k>=4 and key in T[k]]
            vals_all={k: (T[k][key].get(r,(0,0,0,0))[comp]) for k in ks}
            # find first k from which polynomial of degree floor(delta/2) (exactly) fits all later data
            res=None
            for k0 in ks:
                sub=[k for k in ks if k>=k0]
                vals=[Fr(vals_all[k]) for k in sub]
                dg,npts=polydeg(sub,vals)
                if dg is not None and npts>=dg+3:   # at least 2 spare points
                    res=(k0,dg,npts); break
            print(f"  delta={delta} (r={r}) comp={comp}: first k0 with polynomial tail = {res}, floor(delta/2)={delta//2}")
