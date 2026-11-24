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

# construction 1 (step_45 variant A): relations only from (a)sh(bc), summed, NOT deduped
rels1 = []
for a in A:
    for v in itertools.product(A, repeat=2):
        rels1.append(tuple(sh((a,), v)))
# construction 2 (step_46): both orders, deduped
rels2 = set()
for a in A:
    for b,c in itertools.product(A, repeat=2):
        rels2.add(tuple(sorted(tuple(x) for x in sh((a,), (b,c)))))
        rels2.add(tuple(sorted(tuple(x) for x in sh((a,b), (c,)))))
S1 = set(tuple(sorted(r)) for r in rels1)
S2 = set(rels2)
print('|rels1| dedup:', len(S1), ' |rels2|:', len(S2))
print('in S2 not S1:', S2 - S1)
print('in S1 not S2:', S1 - S2)
