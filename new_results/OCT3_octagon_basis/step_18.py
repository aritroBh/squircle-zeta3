import json, time
from mpmath import mp, mpf, sqrt, sin, cos, pi, nstr, nint, zeta
mp.dps = 280
from flint import fmpz_mat

b34 = json.load(open("basis34.json"))
names34 = b34["names"]; X34 = [mpf(s) for s in b34["vals"]]
base17 = X34[::2]; basenames17 = names34[::2]
Pk = [mpf(s) for s in json.load(open("pder.json"))["vals"]]
r2 = sqrt(2); u = sqrt(2+r2); z3 = zeta(3)
mults=[(mpf(1),"1"),(r2,"r2"),(u,"u"),(u*r2,"r2*u")]
X68=[]; names68=[]
for nm,c in zip(basenames17, base17):
    for mv,mn in mults:
        X68.append(c*mv); names68.append(nm if mn=="1" else mn+"*"+nm)

v=[Pk[1]]+X68; vn=["P'(1)"]+names68; n=len(v)
B=10**260
rows=[]
for i in range(n):
    r=[0]*(n+1); r[i]=1; r[n]=int(nint(B*v[i])); rows.append(r)
Rl=[[int(x) for x in row] for row in fmpz_mat(rows).lll().tolist()]
best=None
for row in Rl:
    a=row[:n]
    s=sum(mpf(ai)*xi for ai,xi in zip(a,v) if ai)
    if best is None or abs(s)<best[0]: best=(abs(s),a)
s,a=best
print("residual:", nstr(s,3))
print("relation:")
terms=[(ai,vn[i]) for i,ai in enumerate(a) if ai]
for ai,nm in terms: print("   %+4d * %s" % (ai,nm))
# solve for P'
import sympy as sp
a0=sp.Integer(a[0])
print("\nP'(1) =")
res = -sum(sp.Integer(ai)*sp.Symbol("x%d"%i) for i,ai in enumerate(a[1:],1) if ai)/a0
for i,ai in enumerate(a[1:],1):
    if ai:
        c = sp.Rational(-ai, a0)
        print("   %+s * %s" % (c, vn[i]))
