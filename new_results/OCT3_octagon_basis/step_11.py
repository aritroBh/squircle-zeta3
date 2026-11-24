import json, time
from mpmath import mp, mpf, nstr, nint
mp.dps = 280
from flint import fmpz_mat
import sympy as sp
t0=time.time()

d = json.load(open("consts62.json"))
names = d["names"]; X = [mpf(s) for s in d["vals"]]
n = len(X); B = 10**255
rows = []
for i in range(n):
    row=[0]*(n+1); row[i]=1; row[n]=int(nint(B*X[i])); rows.append(row)
Rl = [[int(v) for v in row] for row in fmpz_mat(rows).lll().tolist()]
rels=[]
for row in Rl:
    a=row[:n]
    s=sum(mpf(ai)*xi for ai,xi in zip(a,X) if ai)
    if abs(s)<mpf('1e-278'): rels.append(a)
A = sp.Matrix(rels)
R, pivots = A.rref()
pivots=list(pivots)
indep = [j for j in range(n) if j not in pivots]
print("rank:", len(pivots), " independent:", len(indep))
print("independent constants:")
for j in indep: print("   ", names[j])

print("\n--- eliminated constants expressed in the basis (each verified) ---")
maxres=mpf(0)
for i,p in enumerate(pivots):
    terms=[(R[i,j],j) for j in indep if R[i,j]!=0]
    # verify: x_p + sum R_ij x_j = 0
    s = X[p] + sum(mpf(str(c))*X[j] for c,j in terms)
    maxres=max(maxres,abs(s))
    rhs = " ".join("%s*%s"%(("%+s"%(-c) if c!=-1 else "+1"), names[j]) for c,j in terms)
    print("%s = %s   |res|=%s" % (names[p], rhs if rhs else "0", nstr(abs(s),1)))
print("max residual over 28 expressions:", nstr(maxres,2))
json.dump({"indep_idx":indep, "names":[names[j] for j in indep],
           "vals":[nstr(X[j],280) for j in indep]}, open("basis34.json","w"))
print("saved basis34.json  [%.1fs]"%(time.time()-t0))
