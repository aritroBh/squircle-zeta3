
import sys
sys.path.insert(0, "<project>/research/08_theory/G2b_point_depth_M4")
from fractions import Fraction as Fr
from g2b_dbg import s05
C0 = Fr(1,2)
w = (0, 1, 0, -1, 1)
S5w = s05.S_J(w, C0, 0)
print("5-symbol of J((0,1,0,-1,1); 0->1/2):")
for k, v in sorted(S5w.items()):
    print("   ", [tuple(str(x) for x in l) for l in k], v)
D5w = s05.delta_minus_5sym(S5w)
print("delta^-(5sym):")
for k, v in D5w.items():
    print("   ", [tuple(str(x) for x in l) for l in k], v)
print("\ndelta_word terms (my port):")
pts = (Fr(0),) + tuple(w) + (C0,)
for k in (2, 3):
    for i in range(6 - k):
        gap = w[i:i + k]; ret = w[:i] + w[i + k:]
        ei, ci = pts[i], pts[i + k + 1]
        if k == 2:
            lam = s05.lam2(gap, ci, ei)
            print(f"  k=2 i={i}: gap={gap} {ei}->{ci} lam={lam}; ret={ret} 0->1/2, S(ret):")
            if lam:
                for lw, cf in s05.S_J(ret, C0, Fr(0)).items():
                    print(f"        {[tuple(str(x) for x in l) for l in lw]} {cf}")
        else:
            lam = s05.lam2(ret, C0, Fr(0))
            print(f"  k=3 i={i}: ret={ret} 0->1/2 lam={lam}; gap={gap} {ei}->{ci}, S(gap):")
            if lam:
                for lw, cf in s05.S_J(gap, ci, ei).items():
                    print(f"        {[tuple(str(x) for x in l) for l in lw]} {cf}")
