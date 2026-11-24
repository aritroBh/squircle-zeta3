import itertools
from flint import fmpz_mat
def sh(u, v):
    if not u: return [v]
    if not v: return [u]
    out = []
    for x in sh(u[1:], v): out.append((u[0],)+x)
    for x in sh(u, v[1:]): out.append((v[0],)+x)
    return out
for s in (1,2,3,4):
    A = ["0","1"]+[f"p{i+1}" for i in range(s)]
    words3 = list(itertools.product(A, repeat=3))
    def good(w): return len(set(w))==3 and w[0] != "0" and w[2] != "1"
    good_words = [w for w in words3 if good(w)]
    bad_words = [w for w in words3 if not good(w)]
    rels = set()
    for a in A:
        for b,c in itertools.product(A, repeat=2):
            rels.add(tuple(sorted(tuple(x) for x in sh((a,),(b,c)))))
            rels.add(tuple(sorted(tuple(x) for x in sh((a,b),(c,)))))
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
    bad_left = []
    for b in bad_words:
        v = [0]*len(words3); v[widx[b]] = 1
        if fmpz_mat(rows+[v]).rank() != rMR: bad_left.append(b)
    conv = [w for w in bad_left if w[0]!="0" and w[2]!="1"]
    print(f"s={s}: NOT reducible: {len(bad_left)} of which CONVERGENT: {len(conv)}")
    if s==1: print('   s=1 not reducible:', bad_left)
    if s==2 and conv: print('   s=2 convergent not reducible:', conv[:20])
