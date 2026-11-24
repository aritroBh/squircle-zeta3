
import sys
sys.path.insert(0, "<project>/research/08_theory/G2b_point_depth_M4")
from fractions import Fraction as Fr
from g2b_dbg import s05
C0 = Fr(1,2)

def delta_word_variant(w, e, c):
    pts = (e,) + tuple(w) + (c,)
    out = {}
    for k in (2, 3):
        for i in range(6 - k):
            gap = w[i:i + k]; ret = w[:i] + w[i + k:]
            ei, ci = pts[i], pts[i + k + 1]
            if k == 2:
                lam = s05.lam2(gap, ei, ci)
                if not lam: continue
                S = s05.S_J(ret, c, e)
                for lw, cf in S.items():
                    out[lw] = out.get(lw, Fr(0)) + cf * lam
            else:
                lam = s05.lam2(ret, c, e)
                if not lam: continue
                S = s05.S_J(gap, ei, ci)
                for lw, cf in S.items():
                    out[lw] = out.get(lw, Fr(0)) - cf * lam
    return {k: v for k, v in out.items() if v}

S5 = s05.omega5_half()
D5 = s05.delta_minus_5sym(S5)
D = {}
for w, coef in s05.signed_words():
    for lw, v in delta_word_variant(w, 0, C0).items():
        D[lw] = D.get(lw, Fr(0)) + coef * v
D = {k: v for k, v in D.items() if v}
keys = sorted(set(D) | set(D5))
ratios = {}
for k in keys:
    a = D.get(k, Fr(0)); b = D5.get(k, Fr(0))
    ratios[k] = (a / b) if b else None
from collections import Counter
print("ratio multiset:", Counter(ratios.values()))
nmiss = [k for k in keys if ratios[k] not in (Fr(1),)]
print("non-ratio-1 count:", len(nmiss), "of", len(keys))
for k in nmiss[:12]:
    print(f"  {k}: variant={D.get(k, Fr(0))}  D5={D5.get(k, Fr(0))}")
