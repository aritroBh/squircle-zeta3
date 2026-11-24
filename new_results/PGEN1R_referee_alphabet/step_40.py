from fractions import Fraction as Fr
import itertools
from flint import fmpz_mat

def sh(u, v):
    if not u: return [v]
    if not v: return [u]
    out = []
    for x in sh(u[1:], v): out.append((u[0],)+x)
    for x in sh(u, v[1:]): out.append((v[0],)+x)
    return out

def check_shuffle_pgen1(s):
    A = ["0","1"]+[f"p{i+1}" for i in range(s)]
    words3 = list(itertools.product(A, repeat=3))
    def good(w):
        return len(set(w))==3 and w[0] not in ("0","1") and w[2] != "1"
    good_words = [w for w in words3 if good(w)]
    bad_words = [w for w in words3 if not good(w)]
    rels = []
    for a in A:
        for b,c in itertools.product(A, repeat=2):
            rels.append([tuple(x) for x in sh((a,),(b,c))])
            rels.append([tuple(x) for x in sh((a,b),(c,))])
    rels = list(set(tuple(sorted(r)) for r in rels))
    widx = {w:i for i,w in enumerate(words3)}
    rows = []
    for r in rels:
        v = [0]*len(words3)
        for w in r: v[widx[w]] = 1
        rows.append(v)
    for g in good_words:
        v = [0]*len(words3); v[widx[g]] = 1
        rows.append(v)
    M = fmpz_mat(rows)
    rMR = M.rank()
    fails = []
    for b in bad_words:
        v = [0]*len(words3); v[widx[b]] = 1
        M2 = fmpz_mat(rows + [v])
        if M2.rank() != rMR:
            fails.append(b)
    print(f"s={s}: #good={len(good_words)}, #bad={len(bad_words)}, #FAIL={len(fails)}")
    if fails and s<=2: print('   failing:', fails[:12])
for s in (1,2,3,4):
    check_shuffle_pgen1(s)
