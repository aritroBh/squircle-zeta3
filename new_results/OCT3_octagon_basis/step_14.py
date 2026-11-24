import json, time
from mpmath import mp, mpf, nstr, nint
mp.dps = 280
from flint import fmpz_mat
t0=time.time()

b34 = json.load(open("basis34.json"))
names = b34["names"]; X = [mpf(s) for s in b34["vals"]]
P3 = mpf(json.load(open("pder.json"))["vals"][3])

v = [P3] + X
vn = ["P'''(1)"] + names
n = len(v)
B = 10**260
rows=[]
for i in range(n):
    r=[0]*(n+1); r[i]=1; r[n]=int(nint(B*v[i])); rows.append(r)
Rl = [[int(x) for x in row] for row in fmpz_mat(rows).lll().tolist()]
print("LLL done [%.1fs]"%(time.time()-t0))
cands=[]
for row in Rl:
    a=row[:n]
    s=sum(mpf(ai)*xi for ai,xi in zip(a,v) if ai)
    cands.append((abs(s), max(abs(ai) for ai in a), abs(a[0]), a))
cands.sort()
print("smallest-residual combos (|s|, maxcoeff, |coeff of P'''|):")
for c in cands[:10]:
    print("  |s|=%s  maxc=%-12d |a0|=%d" % (nstr(c[0],3), c[1], c[2]))
