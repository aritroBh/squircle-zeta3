
import sys
sys.path.insert(0, "<project>/research/08_theory/G2b_point_depth_M4")
from fractions import Fraction as Fr
from g2b_dbg import s05

def delta_word_variant(w, e, c, sign3=-1, gap_rev=False, ret_rev=False, lam_scale=Fr(1)):
    pts = (e,) + tuple(w) + (c,)
    out = {}
    for k in (2, 3):
        for i in range(6 - k):
            gap = w[i:i + k]; ret = w[:i] + w[i + k:]
            ei, ci = pts[i], pts[i + k + 1]
            if k == 2:
                lam = s05.lam2(gap, ei, ci) if gap_rev else s05.lam2(gap, ci, ei)
                lam *= lam_scale
                if not lam: continue
                S = s05.S_J(ret, e, c) if ret_rev else s05.S_J(ret, c, e)
                for lw, cf in S.items():
                    out[lw] = out.get(lw, Fr(0)) + cf * lam
            else:
                lam = s05.lam2(ret, e, c) if ret_rev else s05.lam2(ret, c, e)
                lam *= lam_scale
                if not lam: continue
                S = s05.S_J(gap, ei, ci) if gap_rev else s05.S_J(gap, ci, ei)
                for lw, cf in S.items():
                    out[lw] = out.get(lw, Fr(0)) + sign3 * cf * lam
    return {k: v for k, v in out.items() if v}

w = (0, 1, 0, -1, 1)
C0 = Fr(1,2)
S5w = s05.S_J(w, C0, 0)
D5w = s05.delta_minus_5sym(S5w)
print("D5 has", len(D5w), "words")
for gap_rev in (False, True):
    for ret_rev in (False, True):
        for sign3 in (-1, 1):
            for lam_scale in (Fr(1), Fr(2)):
                D = delta_word_variant(w, 0, C0, sign3, gap_rev, ret_rev, lam_scale)
                if D == D5w:
                    print("MATCH:", gap_rev, ret_rev, sign3, lam_scale)
# print both for the default to inspect
D = delta_word_variant(w, 0, C0)
keys = sorted(set(D) | set(D5w))
for k in keys:
    print(f"  {k}: variant={D.get(k, Fr(0))}  D5={D5w.get(k, Fr(0))}")
