import itertools
from flint import fmpz_mat
# recreate my EARLIER (failing) test verbatim and debug
s = 1
A = ["0","1"]+[f"p{i+1}" for i in range(s)]
def sh(u, v):
    if not u: return [v]
    if not v: return [u]
    out = []
    for x in sh(u[1:], v): out.append((u[0],)+x)
    for x in sh(u, v[1:]): out.append((v[0],)+x)
    return out
words3 = list(itertools.product(A, repeat=3))
widx = {w:i for i,w in enumerate(words3)}
Wdc = [w for w in words3 if len(set(w))==3 and w[0]!="0" and w[2]!="1"]
print('Wdc:', Wdc)
# variant A: relations as SUM vectors (my earlier test)
rowsA = []
for w in Wdc:
    r=[0]*27; r[widx[w]]=1; rowsA.append(r)
for a in A:
    for v in itertools.product(A, repeat=2):
        r = [0]*27
        for w in sh((a,), v):
            r[widx[w]] += 1
        rowsA.append(r)
rA = fmpz_mat(rowsA).rank()
# variant B: PGEN1 — note they add BOTH (a)sh(bc) and (ab)sh(c), dedup by sorted
rels = set()
for a in A:
    for b,c in itertools.product(A, repeat=2):
        rels.add(tuple(sorted(tuple(x) for x in sh((a,),(b,c)))))
        rels.add(tuple(sorted(tuple(x) for x in sh((a,b),(c,))))))
rowsB = []
for g in Wdc:
    r=[0]*27; r[widx[g]]=1; rowsB.append(r)
for rr in rels:
    r=[0]*27
    for w in rr: r[widx[w]] += 1
    rowsB.append(r)
rB = fmpz_mat(rowsB).rank()
print('rank variant A (only (a)sh(bc), summed):', rA)
print('rank variant B (both, dedup-sorted):', rB)
# check (1,1,p1) membership both ways
w0 = ('1','1','p1')
v0 = [0]*27; v0[widx[w0]]=1
print('(1,1,p1) in span A:', fmpz_mat(rowsA+[v0]).rank()==rA, ' in span B:', fmpz_mat(rowsB+[v0]).rank()==rB)
