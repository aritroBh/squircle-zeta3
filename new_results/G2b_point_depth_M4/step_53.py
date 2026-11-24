
import sys
sys.path.insert(0, "<project>/research/08_theory/G2b_point_depth_M4")
from fractions import Fraction as Fr
from g2b_dbg import s05
C0 = Fr(1,2)

def delta_word_variant(w, e, c, sign3, gap_rev, ret_rev, lam_scale):
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

S5 = s05.omega5_half()
D5 = s05.delta_minus_5sym(S5)
print("D5 words:", len(D5))
words = s05.signed_words()
for gap_rev in (False, True):
    for ret_rev in (False, True):
        for sign3 in (-1, 1):
            D = {}
            for w, coef in words:
                for lw, v in delta_word_variant(w, 0, C0, sign3, gap_rev, ret_rev, Fr(1)).items():
                    D[lw] = D.get(lw, Fr(0)) + coef * v
            D = {k: v for k, v in D.items() if v}
            same = D == D5
            neg = all(D.get(k, Fr(0)) == -v for k, v in D5.items()) and all(-v == D5.get(k) for k, v in D.items())
            print(f"gap_rev={gap_rev} ret_rev={ret_rev} sign3={sign3}: words={len(D)} match={same} minus-match={neg}")
