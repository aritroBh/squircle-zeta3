import os, sys, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
import point_engine as PE, mzv_exact as M, kimilib
from fractions import Fraction as F
pt=(2,3,5,7,11)
t0=time.time(); vw=PE.volreg_words_point(pt,12); t1=time.time()
print(f"volreg_words_point: {t1-t0:.2f}s, words at W12: {len(vw[12])}")
t0=time.time(); m0=kimilib.m0_coords(pt,12); t1=time.time()
print(f"m0_coords: {t1-t0:.2f}s")
print("m0 W12 coords 7..11:", m0[12][7:])
print("m0 CLASS:", m0[12][9]+F(14,9)*m0[12][10])
# full reduced vs vw+vwr only, coords 7..11
c=PE.reduced_at(pt,12)[12]
T=PE.table(12)
vwr=PE.volreg_words_point(pt[::-1],12)
acc=[F(0)]*12
for src in (vw[12],vwr[12]):
    for w,pp in src.items():
        rw=T.reduce_word(w)
        for b in range(12): acc[b]+=pp*rw[b]
print("full reduced 7..11:   ", c[7:])
print("CLASS full:", c[9]+F(14,9)*c[12-2])  # c[10]
print("vw+vwr only 7..11:    ", acc[7:])
print("CLASS vw+vwr:", acc[9]+F(14,9)*acc[10])
