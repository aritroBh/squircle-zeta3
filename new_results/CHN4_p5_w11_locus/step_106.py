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
def w12irr(pt):
    c12=coords_reg_idx(pt,12)
    CLASS=c12[9]+F(14,9)*c12[10]
    return (c12[7],c12[8],CLASS,c12[11])
print("W12 irreducibles (z2^2z35, z2z37, CLASS, z2118) on faces/strata:")
for lbl,pts in [
   ('a2=0',[(1,0,2,3,1),(2,0,1,1,3),(3,0,2,1,2)]),
   ('a3=0',[(1,2,0,3,1),(2,1,0,1,3)]),
   ('a4=0',[(1,2,3,0,1),(2,1,1,0,3)]),
   ('a5=0 stratum {a2=a4,a1=a3}',[(1,2,1,2,0),(2,3,2,3,0),(1,1,1,1,0)]),
]:
    for pt in pts:
        r=w12irr(pt)
        print(f"  {lbl} {pt}: {r}  {'ALL ZERO' if all(x==0 for x in r) else 'NONZERO'}")
