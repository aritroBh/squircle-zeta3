import sys, os
sys.dont_write_bytecode=True
sys.path.insert(0, os.getcwd())
from chk_feed import *
from math import comb
K={}
nz_i=0; nz_c=0
for k in range(6,53):
    W=kernels_deg5(towers_at(k))
    K[k]={}
    for key in KEYS5:
        simple,higher,const=pf_split(W[key])
        if any(not x.is_zero() for x in simple[1]) or any(not x.is_zero() for x in simple[2]): nz_i+=1
        if not const.is_zero(): nz_c+=1
        K[k][key]={p:c for p,c in realdict(simple).items()}
print('kernels built from degree-6 towers for k=6..52; nonzero +-i simple poles:',nz_i,' nonzero const:',nz_c)
def vec(k,key,p): return K[k][key].get(p,(Fr(0),)*4)
poles=set()
for k in K:
    for key in KEYS5: poles|=set(K[k][key].keys())
print('pole union', sorted((str(a),str(b)) for a,b in poles))
tot=bad=0
for key in KEYS5:
    for p in poles:
        for k0 in range(6,48):
            for t in range(4):
                s=sum((-1)**(5-j)*comb(5,j)*vec(k0+j,key,p)[t] for j in range(6))
                tot+=1; bad+=(s!=0)
print('5th differences over every window k0..k0+5, k0=6..47 (so k up to 52), 3 keys x poles x 4 comps:',tot,'bad',bad)
# degree exactly 4 on the whole range 6..52 for every key
for key in KEYS5:
    d4=False
    for p in poles:
        for t in range(4):
            s=sum((-1)**(4-j)*comb(4,j)*vec(6+j,key,p)[t] for j in range(5))
            if s!=0: d4=True
    print(key,'4th difference nonzero (degree exactly 4):',d4)
