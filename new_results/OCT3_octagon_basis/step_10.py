import json, time
from mpmath import mp, mpf, nstr, nint
mp.dps = 280
from flint import fmpz_mat
import sympy as sp

d = json.load(open("consts62.json"))
names = d["names"]; X = [mpf(s) for s in d["vals"]]
n = len(X); B = 10**255
rows = []
for i in range(n):
    row = [0]*(n+1); row[i]=1; row[n] = int(nint(B*X[i])); rows.append(row)
Rl = [[int(v) for v in row] for row in fmpz_mat(rows).lll().tolist()]
rels = []
for row in Rl:
    a = row[:n]
    s = sum(mpf(ai)*xi for ai,xi in zip(a,X) if ai)
    if abs(s) < mpf('1e-278'):
        rels.append((a,s))
print("true relations found:", len(rels))

def rel_str(a):
    terms=[]
    for ai,nm in zip(a,names):
        if ai: terms.append(("%+d"%ai)+"*"+nm)
    return " ".join(terms)
for a,s in rels:
    print("RES", nstr(s,2), "|", rel_str(a))

# rref over QQ
A = sp.Matrix([a for a,_ in rels])
R, pivots = A.rref()
print("\nrref pivots:", pivots, " rank:", len(pivots))
