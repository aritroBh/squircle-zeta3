import json, time
from mpmath import mp, mpf, sqrt, nstr, nint
mp.dps = 280
from flint import fmpz_mat
t0=time.time()

b34 = json.load(open("basis34.json"))
names34 = b34["names"]; X34 = [mpf(s) for s in b34["vals"]]
base17 = X34[::2]; basenames17 = names34[::2]
r2 = sqrt(2); u = sqrt(2+r2)
print("u = sqrt(2+sqrt2) =", nstr(u, 30))
print("C4 check: (54-39*sqrt2)*u =", nstr((54-39*r2)*u, 25))
C4 = 192*( (2-r2)/4 )**0 # placeholder
from mpmath import sin, cos, pi
C4 = 192*sin(pi/8)**4*cos(pi/8)*(3*sin(pi/8)**2-1)
print("C4 =", nstr(C4, 25), " diff:", nstr(C4-(54-39*r2)*u, 5))

mults = [(mpf(1),"1"), (r2,"r2"), (u,"u"), (u*r2,"r2*u")]
X68=[]; names68=[]
for nm,c in zip(basenames17, base17):
    for mv,mn in mults:
        X68.append(c*mv); names68.append(nm if mn=="1" else mn+"*"+nm)
print("68 constants built")

def lll_rows(vals, B):
    n=len(vals); rows=[]
    for i in range(n):
        r=[0]*(n+1); r[i]=1; r[n]=int(nint(B*vals[i])); rows.append(r)
    return [[int(x) for x in row] for row in fmpz_mat(rows).lll().tolist()]

# internal relations among the 68
Rl = lll_rows(X68, 10**260)
cs=[]
for row in Rl:
    a=row[:68]
    s=sum(mpf(ai)*xi for ai,xi in zip(a,X68) if ai)
    cs.append((abs(s), max(abs(ai) for ai in a)))
cs.sort()
print("internal check on 68: smallest |s|:")
for c in cs[:5]: print("  |s|=%s maxc=%d" % (nstr(c[0],3), c[1]))

# hunt P''' against the 68
P3 = mpf(json.load(open("pder.json"))["vals"][3])
v = [P3]+X68; vn=["P'''(1)"]+names68
Rl = lll_rows(v, 10**260)
cs=[]
for row in Rl:
    a=row[:69]
    s=sum(mpf(ai)*xi for ai,xi in zip(a,v) if ai)
    cs.append((abs(s), max(abs(ai) for ai in a), abs(a[0]), a))
cs.sort()
print("\nP''' hunt vs 68 constants: smallest |s|:")
for c in cs[:10]: print("  |s|=%s maxc=%-10d |a0|=%d" % (nstr(c[0],3), c[1], c[2]))
print("[%.1fs]"%(time.time()-t0))
