import itertools
from flint import fmpz_mat
def sh(u, v):
    if not u: return [v]
    if not v: return [u]
    out = []
    for x in sh(u[1:], v): out.append((u[0],)+x)
    for x in sh(u, v[1:]): out.append((v[0],)+x)
    return out
s = 1
A = ["0","1"]+[f"p{i+1}" for i in range(s)]
words3 = list(itertools.product(A, repeat=3))
widx = {w:i for i,w in enumerate(words3)}
Wdc = [w for w in words3 if len(set(w))==3 and w[0]!="0" and w[2]!="1"]
rels2 = set()
for a in A:
    for b,c in itertools.product(A, repeat=2):
        rels2.add(tuple(sorted(tuple(x) for x in sh((a,), (b,c)))))
        rels2.add(tuple(sorted(tuple(x) for x in sh((a,b), (c,)))))
rowsU = []
for g in Wdc:
    r=[0]*27; r[widx[g]]=1; rowsU.append(r)
rowsR = []
for rr in rels2:
    r=[0]*27
    for w in rr: r[widx[w]] += 1
    rowsR.append(r)
print('rank units alone:', fmpz_mat(rowsU).rank())
print('rank relations alone:', fmpz_mat(rowsR).rank())
print('rank units+relations:', fmpz_mat(rowsU+rowsR).rank())
# inspect the relations: maybe duplicate words within a relation are collapsed by set-sorted dedup?
for rr in sorted(rels2)[:8]:
    print(rr)
