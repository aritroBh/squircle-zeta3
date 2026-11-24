import json, time
from mpmath import mp, mpf, nstr, nint
mp.dps = 280
t0=time.time()
from flint import fmpz_mat

d = json.load(open("consts62.json"))
names = d["names"]; X = [mpf(s) for s in d["vals"]]
n = len(X)
B = 10**255
rows = []
for i in range(n):
    row = [0]*(n+1); row[i]=1
    row[n] = int(nint(B*X[i]))
    rows.append(row)
print("rounding check:", nstr(abs(X[0] - mpf(rows[0][n])/B), 3))
R = fmpz_mat(rows).lll()
Rl = [[int(v) for v in row] for row in R.tolist()]
print("LLL done [%.1fs]"%(time.time()-t0))

cands=[]
for row in Rl:
    a = row[:n]; Mlast = row[n]
    s = mpf(0)
    for ai,xi in zip(a,X):
        if ai: s += mpf(ai)*xi
    cands.append((abs(s), max(abs(ai) for ai in a), a, Mlast))
cands.sort(key=lambda t: t[0])
print("smallest residuals among 62 reduced rows:")
for c in cands[:32]:
    print("  |s|=%s  maxc=%-6d M=%d" % (nstr(c[0],3), c[1], c[3]))
print("  --- gap ---")
for c in cands[32:36]:
    print("  |s|=%s  maxc=%-6d M=%d" % (nstr(c[0],3), c[1], c[3]))
