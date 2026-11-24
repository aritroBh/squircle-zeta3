import json, time
from mpmath import mp, mpf, nstr
mp.dps = 280
t0=time.time()

# test python-flint LLL API
from flint import fmpz_mat
M = fmpz_mat([[1, 10**30],[0, 10**30]])
try:
    R = M.lll()
    print("lll returns:", type(R), R.tolist())
except Exception as e:
    print("lll() failed:", e)

d = json.load(open("consts62.json"))
names = d["names"]; X = [mpf(s) for s in d["vals"]]
n = len(X)
B = 10**255
rows = []
for i in range(n):
    row = [0]*n + [int((B*X[i]).round())] if hasattr((B*X[i]),'round') else None
    # safer rounding:
    row = [0]*(n+1); row[i]=1
    row[n] = int(mpf(B*X[i]).__round__())
    rows.append(row)
M = fmpz_mat(rows)
print("built %dx%d matrix  [%.1fs]"%(n,n+1,time.time()-t0))
R = M.lll()
print("LLL done  [%.1fs]"%(time.time()-t0))
Rl = R.tolist()
print("reduced rows:", len(Rl))
# analyze
cands=[]
for row in Rl:
    a = row[:n]; Mlast = row[n]
    s = sum(mpf(ai)*xi for ai,xi in zip(a,X) if ai)
    cands.append((abs(s), max(abs(ai) for ai in a), a, Mlast, s))
cands.sort(key=lambda t: t[0])
print("smallest |residual| rows:")
for c in cands[:20]:
    print("  |s|=%s  maxc=%d  M=%d" % (nstr(c[0],3), c[1], c[3]))
