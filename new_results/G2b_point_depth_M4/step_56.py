
import sys
sys.path.insert(0, "<project>/research/08_theory/G2b_point_depth_M4")
from fractions import Fraction as Fr
from g2b_dbg import s05
C0 = Fr(1,2)
import g2_fit as gf
from step_10 import load_terms

words5 = []
for kind, w, e, q in load_terms(['hpl_fit_k4_(0,0,0).json']):
    if len(w) != 5: continue
    words5.append(w)
words5 = sorted(set(words5))
print(len(words5), "distinct unsigned words")

stats = {"equal": 0, "minus": 0, "both_zero": 0, "other": 0}
for w in words5:
    D5w = s05.delta_minus_5sym(s05.S_J(w, C0, 0))
    Dw = s05.delta_word_sym(w, 0, C0)
    if not D5w and not Dw:
        stats["both_zero"] += 1
    elif D5w == Dw:
        stats["equal"] += 1
    elif all(Dw.get(k, Fr(0)) == -v for k, v in D5w.items()) and all(-v == D5w.get(k) for k, v in Dw.items()):
        stats["minus"] += 1
    else:
        stats["other"] += 1
        if stats["other"] <= 3:
            print("OTHER:", w)
            for k in sorted(set(Dw) | set(D5w)):
                print("   ", k, Dw.get(k, Fr(0)), D5w.get(k, Fr(0)))
print(stats)
