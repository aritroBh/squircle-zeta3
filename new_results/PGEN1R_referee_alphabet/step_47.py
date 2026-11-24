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
def good(w): return len(set(w))==3 and w[0] != "0" and w[2] != "1"
good_words = [w for w in words3 if good(w)]
rels = set()
for a in A:
    for b,c in itertools.product(A, repeat=2):
        rels.add(tuple(sorted(tuple(x) for x in sh((a,), (b,c)))))
        rels.add(tuple(sorted(tuple(x) for x in sh((a,b), (c,)))))
widx = {w:i for i,w in enumerate(words3)}
rows = []
for r in rels:
    v = [0]*len(words3)
    for w in r: v[widx[w]] = 1
    rows.append(v)
for g in good_words:
    v = [0]*len(words3); v[widx[g]] = 1
    rows.append(v)
rMR = fmpz_mat(rows).rank()
print('rMR =', rMR)
w0 = ('1','1','p1'); v0=[0]*27; v0[widx[w0]]=1
print('(1,1,p1) reducible:', fmpz_mat(rows+[v0]).rank()==rMR)
# inspect a relation containing (1,1,p1)
for rr in rels:
    if w0 in rr: print('rel:', rr)
