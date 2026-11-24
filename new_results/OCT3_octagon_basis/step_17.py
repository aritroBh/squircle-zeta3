import json, time
from mpmath import mp, mpf, sqrt, sin, cos, pi, nstr, nint, zeta
mp.dps = 280
from flint import fmpz_mat
t0=time.time()

b34 = json.load(open("basis34.json"))
names34 = b34["names"]; X34 = [mpf(s) for s in b34["vals"]]
base17 = X34[::2]; basenames17 = names34[::2]
Pk = [mpf(s) for s in json.load(open("pder.json"))["vals"]]
r2 = sqrt(2); u = sqrt(2+r2); z3 = zeta(3)
C4 = 192*sin(pi/8)**4*cos(pi/8)*(3*sin(pi/8)**2-1)
mults=[(mpf(1),"1"),(r2,"r2"),(u,"u"),(u*r2,"r2*u")]
X68=[]; names68=[]
for nm,c in zip(basenames17, base17):
    for mv,mn in mults:
        X68.append(c*mv); names68.append(nm if mn=="1" else mn+"*"+nm)

def hunt(target, tname, vals, vnames, B=10**260, top=6):
    v=[target]+vals; n=len(v)
    rows=[]
    for i in range(n):
        r=[0]*(n+1); r[i]=1; r[n]=int(nint(B*v[i])); rows.append(r)
    Rl=[[int(x) for x in row] for row in fmpz_mat(rows).lll().tolist()]
    cs=[]
    for row in Rl:
        a=row[:n]
        s=sum(mpf(ai)*xi for ai,xi in zip(a,v) if ai)
        cs.append((abs(s), max(abs(ai) for ai in a), abs(a[0]), a))
    cs.sort()
    print("hunt %s (dim %d):" % (tname, n))
    for c in cs[:top]:
        print("   |s|=%s maxc=%-10d |a0|=%d" % (nstr(c[0],3), c[1], c[2]))
    return cs

cs2 = hunt(Pk[2], "P''(1) vs 68", X68, names68)
print()
cs1 = hunt(Pk[1], "P'(1) vs 68", X68, names68)
print()
T = Pk[3] - C4*z3
csT = hunt(T, "P'''-C4*z3 vs 68", X68, names68)
print("[%.1fs]"%(time.time()-t0))
