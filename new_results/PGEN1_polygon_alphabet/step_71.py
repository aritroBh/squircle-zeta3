# step_43: corrected good words: s=4 rerun + final shuffle-reduction check (all words -> good + products)
exec(open("symbolproof_engine.py").read())
import time
t0=time.time()
run_proof(4)
print("s=4 time:", round(time.time()-t0,1),"s")

# shuffle check with corrected good set: good = distinct, a1 != 0, a3 != 1
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
            rels.add(tuple(sorted(tuple(x) for x in sh((a,b),(c,))))))
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
    print(f"s={s}: #words={len(words3)} good={len(good_words)} bad={len(bad_words)}; not reducible: {bad_left}")
