import sys
sys.path.insert(0,'.')
import point_engine as PE, mzv_exact as M
from fractions import Fraction as F
def coords_reg_idx(a, W):
    a=[F(x) for x in a]
    vw=PE.volreg_words_point(a,W); vwr=PE.volreg_words_point(a[::-1],W)
    T=PE.table(W); d=T.nullity; acc=[F(0)]*d
    for src in (vw[W],vwr[W]):
        for w,pp in src.items():
            rw=T.reduce_word(w)
            for b in range(d):
                if rw[b]: acc[b]+=pp*rw[b]
    return acc
for pt,lbl in [((0,1,1,1,1),'t=1,a3=1'), ((0,2,1,2,1),'t=2,a3=1'), ((0,3,2,3,2),'t=3,a3=2')]:
    print(f"\n=== boundary pt {pt} ({lbl}) ===")
    c11=coords_reg_idx(pt,11); c12=coords_reg_idx(pt,12)
    CLASS=c12[9]+F(14,9)*c12[10]
    print(f"  W11: z3z35(idx7)={c11[7]}, f11=z(3,3,5)(idx8)={c11[8]}")
    print(f"  W12: z2^2z35={c12[7]}, z2z37={c12[8]}, z(3,9)={c12[9]}, z(5,7)={c12[10]}, z2118={c12[11]}")
    print(f"  CLASS(z39+14/9 z57)={CLASS}")
