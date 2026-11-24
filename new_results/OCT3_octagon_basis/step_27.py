import json, time
from mpmath import mp, mpf, sqrt, nstr, nint
mp.dps = 520
from flint import fmpz_mat
t0=time.time()
d=json.load(open("hidata.json"))
Pk=[mpf(s) for s in d["vals"]]
C={k:mpf(v) for k,v in d["consts"].items()}
r2=sqrt(2); u=sqrt(2+r2)
w3n=["Li3(2b)","Li3((2+r2)/4)","Li3(12r2-16)","z3","pi2*l2","pi2*lu","l2^3","l2^2*lu","l2*lu^2","lu^3"]
w2n=["pi2","l2^2","l2*lu","lu^2","Li2(12r2-16)"]
w10n=["l2","lu","1"]
basen=w3n+w2n+w10n
XB=[];NB=[]
for nm in basen:
    c=C[nm]
    XB += [c, r2*c]; NB += [nm, "r2*"+nm]
T = Pk[3]/u
v=[T]+XB; n=len(v); B=10**480
rows=[]
for i in range(n):
    r=[0]*(n+1); r[i]=1; r[n]=int(nint(B*v[i])); rows.append(r)
Rl=[[int(x) for x in row] for row in fmpz_mat(rows).lll().tolist()]
print("LLL done [%.1fs]"%(time.time()-t0))
cs=[]
for row in Rl:
    a=row[:n]
    s=sum(mpf(ai)*xi for ai,xi in zip(a,v) if ai)
    cs.append((abs(s),max(abs(ai) for ai in a),abs(a[0]),a))
cs.sort(key=lambda t:t[0])
print("divided hunt P'''/u vs 36, dim 37, B=1e480:")
for c in cs[:8]: print("   |s|=%s maxc=%-14d |a0|=%d"%(nstr(c[0],3),c[1],c[2]))
