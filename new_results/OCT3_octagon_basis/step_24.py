import json, time
from mpmath import mp, mpf, sqrt, log, pi, zeta, polylog, nstr, nint
mp.dps = 280
from flint import fmpz_mat
import sympy as sp
r2=sqrt(2); b=r2-1; u=sqrt(2+r2)
l2=log(2); lu=log(1+r2); pi2=pi**2
Pk=[mpf(s) for s in json.load(open("pder.json"))["vals"]]
Li2g=polylog(2,12*r2-16)
w2b=[pi2,l2**2,l2*lu,lu**2,Li2g,l2,lu,mpf(1)]
w2bn=["pi2","l2^2","l2*lu","lu^2","Li2(12r2-16)","l2","lu","1"]
mults=[(mpf(1),"1"),(r2,"r2"),(u,"u"),(u*r2,"r2*u")]
XB=[];NB=[]
for nm,c in zip(w2bn,w2b):
    for mv,mn in mults: XB.append(c*mv); NB.append(nm if mn=="1" else mn+"*"+nm)
v=[Pk[2]]+XB;n=len(v);B=10**260
rows=[]
for i in range(n):
    r=[0]*(n+1);r[i]=1;r[n]=int(nint(B*v[i]));rows.append(r)
Rl=[[int(x) for x in row] for row in fmpz_mat(rows).lll().tolist()]
best=None
for row in Rl:
    a=row[:n]
    s=sum(mpf(ai)*xi for ai,xi in zip(a,v) if ai)
    if best is None or abs(s)<abs(best[0]): best=(s,a)
s,a=best
print("residual:",nstr(s,3)," maxc:",max(abs(ai) for ai in a))
vn=["P''"]+NB
print("relation terms:")
a0=sp.Integer(a[0])
print("  a0 (P'' coeff) =", a0)
for i,ai in enumerate(a):
    if ai: print("   %+6d * %s"%(ai,vn[i]))
print("\nP''(1) =")
for i,ai in enumerate(a[1:],1):
    if ai: print("   %+s * %s"%(sp.Rational(-ai,a0),vn[i]))
# verify numerically
val = sum(sp.Rational(-ai,a0)*mpf(str(XB[i-1])) for i,ai in enumerate(a[1:],1) if ai)
print("check:", nstr(val-Pk[2],5))
