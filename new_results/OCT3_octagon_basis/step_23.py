import json, time
from mpmath import mp, mpf, sqrt, log, pi, zeta, polylog, nstr, nint
mp.dps = 280
from flint import fmpz_mat
import sympy as sp
t0=time.time()
r2=sqrt(2); b=r2-1; u=sqrt(2+r2)
l2=log(2); lu=log(1+r2); z3=zeta(3); pi2=pi**2
Pk=[mpf(s) for s in json.load(open("pder.json"))["vals"]]
Li2g = polylog(2, 12*r2-16)

def hunt(target,tname,vals,vnames,B=10**260,top=4):
    v=[target]+vals;n=len(v)
    rows=[]
    for i in range(n):
        r=[0]*(n+1);r[i]=1;r[n]=int(nint(B*v[i]));rows.append(r)
    Rl=[[int(x) for x in row] for row in fmpz_mat(rows).lll().tolist()]
    cs=[]
    for row in Rl:
        a=row[:n]
        s=sum(mpf(ai)*xi for ai,xi in zip(a,v) if ai)
        cs.append((abs(s),max(abs(ai) for ai in a),abs(a[0]),a))
    cs.sort(key=lambda t:t[0])
    print("hunt %s (dim %d):"%(tname,n))
    for c in cs[:top]: print("   |s|=%s maxc=%-10d |a0|=%d"%(nstr(c[0],3),c[1],c[2]))
    return cs

# pure weight-2 basis x {1,r2,u,r2u}
w2b=[pi2,l2**2,l2*lu,lu**2,Li2g]
w2bn=["pi2","l2^2","l2*lu","lu^2","Li2(12r2-16)"]
mults=[(mpf(1),"1"),(r2,"r2"),(u,"u"),(u*r2,"r2*u")]
XB=[];NB=[]
for nm,c in zip(w2bn,w2b):
    for mv,mn in mults: XB.append(c*mv); NB.append(nm if mn=="1" else mn+"*"+nm)
cs2=hunt(Pk[2],"P'' pure-w2",XB,NB)
if cs2[0][0]<mpf('1e-250'):
    a=cs2[0][3]; s=cs2[0][0]
    print("\nP'' relation (residual %s):"%nstr(s,3))
    a0=sp.Integer(a[0])
    for i,ai in enumerate(a):
        if ai: print("   %+6d * %s"%(ai, (["P''"]+NB)[i]))
    print("\nP''(1) =")
    for i,ai in enumerate(a[1:],1):
        if ai: print("   %+s * %s"%(sp.Rational(-ai,a0), NB[i]))
print("[%.1fs]"%(time.time()-t0))
