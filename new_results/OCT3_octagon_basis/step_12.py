import json, time
from mpmath import mp, mpf, nstr, nint, pslq
mp.dps = 280
from flint import fmpz_mat
t0=time.time()

d = json.load(open("basis34.json"))
names = d["names"]; X = [mpf(s) for s in d["vals"]]
n = len(X)
B = 10**260
rows=[]
for i in range(n):
    r=[0]*(n+1); r[i]=1; r[n]=int(nint(B*X[i])); rows.append(r)
Rl = [[int(v) for v in row] for row in fmpz_mat(rows).lll().tolist()]
cands=[]
for row in Rl:
    a=row[:n]
    s=sum(mpf(ai)*xi for ai,xi in zip(a,X) if ai)
    cands.append((abs(s), max(abs(ai) for ai in a)))
cands.sort()
print("completeness: smallest-residual integer combos on the 34 basis constants:")
for c in cands[:6]:
    print("  |s|=%s  maxcoeff=%d" % (nstr(c[0],3), c[1]))
print("  [%.1fs]"%(time.time()-t0))

print("\nrunning mpmath.pslq (tol=1e-250, maxcoeff=1e12) ...")
r = pslq(X, tol=mpf('1e-250'), maxcoeff=10**12, maxsteps=100000)
print("pslq result:", r, " [%.1fs]"%(time.time()-t0))
