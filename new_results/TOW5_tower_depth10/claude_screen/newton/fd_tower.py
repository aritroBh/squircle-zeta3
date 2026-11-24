import pickle, random, sys
from fractions import Fraction as Fr
W="<project>/research/08_theory/TOW5_tower_depth10"
which=sys.argv[1] if len(sys.argv)>1 else "K24"
fn={"K24":"data_towers_D10_K24.pkl","K30":"data_towers_D10_K30_deg18.pkl"}[which]
T=pickle.load(open(W+"/"+fn,"rb"))
ks=sorted(T)
print("levels",ks[0],"..",ks[-1], "n levels",len(ks))
def diffs(seq,order):
    s=list(seq)
    for _ in range(order): s=[b-a for a,b in zip(s,s[1:])]
    return s
random.seed(1)
cells=[]  # (key,delta)
keyset=set()
for k in ks:
    for key,dd in T[k].items():
        for dl in dd: keyset.add((key,dl))
cells=sorted(keyset)
print("total (key,delta) series:",len(cells))
n_ok=n_bad=0; exact_deg=0; tested_pts=0; short=0; lower=0
bad=[]
for (key,dl) in cells:
    p,q=key; n=p+q; diag=(p==q)
    k0=n+1 if diag else n
    if dl>2*n-4: continue
    ser=[(k,T[k][key][dl]) for k in ks if k>=k0 and key in T[k] and dl in T[k][key]]
    kk=[k for k,_ in ser]
    if kk!=list(range(k0,k0+len(kk))): bad.append((key,dl,"nonconsecutive")); continue
    P=dl//2
    if len(ser)<P+3: short+=1; continue
    ncomp=len(ser[0][1])
    top_nonzero=False
    for c in range(ncomp):
        s=[Fr(v[c]) for _,v in ser]
        if any(x!=0 for x in diffs(s,P+1)): bad.append((key,dl,c)); break
        else:
            if any(x!=0 for x in diffs(s,P)): top_nonzero=True
    else:
        n_ok+=1; tested_pts+=len(ser)
        if top_nonzero: exact_deg+=1
        else: lower+=1
print("series OK (degree <= floor(delta/2), all components):",n_ok," with exact degree floor(delta/2):",exact_deg," lower degree:",lower," too short:",short)
print("BAD:",len(bad),bad[:10])
# negative control: perturb one value in a random series and confirm detection
cnt=0;det=0
pop=[c for c in cells if c[1]>=2 and c[1]<=2*(c[0][0]+c[0][1])-4]
for (key,dl) in random.sample(pop,min(200,len(pop))):
    p,q=key;n=p+q;k0=n+1 if p==q else n
    ser=[Fr(T[k][key][dl][3]) for k in ks if k>=k0 and key in T[k] and dl in T[k][key]]
    P=dl//2
    if len(ser)<P+3: continue
    j=random.randrange(len(ser)); ser[j]+=Fr(1,10**6)
    cnt+=1
    if any(x!=0 for x in diffs(ser,P+1)): det+=1
print("negative control: perturbed",cnt,"series by 1e-6 in one entry; detected",det)
# also test k = k0-1 for off-diagonal? (should typically fail = initial level exceptional only for diag k=2p)
