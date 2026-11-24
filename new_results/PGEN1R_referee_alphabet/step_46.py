import itertools
from flint import fmpz_mat
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
rels = set()
for a in A:
    for b,c in itertools.product(A, repeat=2):
        rels.add(tuple(sorted(sh((a,), (b,c)))))
        rels.add(tuple(sorted(sh((a,b), (c,)))))
rowsB = []
for g in Wdc:
    r=[0]*27; r[widx[g]]=1; rowsB.append(r)
for rr in rels:
    r=[0]*27
    for w in rr: r[widx[w]] += 1
    rowsB.append(r)
rB = fmpz_mat(rowsB).rank()
print('rank A (a)sh(bc) only:', rA, ' rank B (both, dedup):', rB)
w0 = ('1','1','p1'); v0=[0]*27; v0[widx[w0]]=1
print('(1,1,p1) in A:', fmpz_mat(rowsA+[v0]).rank()==rA, ' in B:', fmpz_mat(rowsB+[v0]).rank()==rB)
print('sh((1,),(p,0)) =', sh(('1',),('p','0')))
print('sh((1,p),(0,)) =', sh(('1','p'),('0',)))
