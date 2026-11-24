import time, sys, pickle
sys.path.insert(0,'.')
import point_engine as PE, mzv_exact as M, kimilib
from fractions import Fraction as F

def coords_reg(a, W):
    """regular part only (vw + reversed), reduced to basis, as Fractions."""
    vw = PE.volreg_words_point(a, W)
    vwr = PE.volreg_words_point(a[::-1], W)
    T = PE.table(W); d=T.nullity
    acc=[F(0)]*d
    for src in (vw[W], vwr[W]):
        for w,pp in src.items():
            rw=T.reduce_word(w)
            for b in range(d):
                if rw[b]: acc[b]+=pp*rw[b]
    return acc

def coords_full(a, W):
    """regular + m0."""
    reg = coords_reg(a,W)
    m0 = kimilib.m0_coords(a,W)[W]
    return [reg[b]+m0[b] for b in range(len(reg))]

pt=(1,2,3,4,5)
ref = pickle.load(open('test_pt_w1112.pkl','rb'))
print("ref pt:", ref['pt'])

for W in (11,):
    reg = coords_reg(pt,W)
    full = coords_full(pt,W)
    print(f"\nW={W} basis:")
    for i,mon in enumerate(M.BASIS[W]):
        lbl='*'.join(','.join(map(str,w)) for w in mon)
        print(f"  [{i}] {lbl:14s} reg={reg[i]}  full={full[i]}  ref={ref['c11'][i]}  match_full={full[i]==ref['c11'][i]}")
