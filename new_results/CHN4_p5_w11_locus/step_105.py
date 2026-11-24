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
tests={
 'a2=0':(1,0,1,1,1), 'a3=0':(1,1,0,1,1), 'a4=0':(1,1,1,0,1),
 'a1=0 off-stratum':(0,1,1,2,1), 'a5=0 off-stratum':(1,1,2,1,0),
}
for lbl,pt in tests.items():
    c12=coords_reg_idx(pt,12); c11=coords_reg_idx(pt,11); c8=coords_reg_idx(pt,8); c10=coords_reg_idx(pt,10)
    CLASS=c12[9]+F(14,9)*c12[10]
    print(f"{lbl} {pt}:")
    print(f"   W8z35={c8[3]} W10z37={c10[6]} f11={c11[8]}")
    print(f"   W12: z2^2z35={c12[7]} z2z37={c12[8]} CLASS={CLASS} z2118={c12[11]}")
