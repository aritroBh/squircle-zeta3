import time, sys, pickle
sys.path.insert(0,'.')
import point_engine as PE, mzv_exact as M, kimilib
from fractions import Fraction as F

def coords_reg(a, W):
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

# Verify W12 at (1,2,3,4,5) against reference (takes a bit longer)
pt=(1,2,3,4,5)
ref = pickle.load(open('test_pt_w1112.pkl','rb'))
t0=time.time()
reg12 = coords_reg(pt,12)
t1=time.time()
print(f"W12 reg computed in {t1-t0:.1f}s")
for i in (7,8,9,10,11):
    print(f"  [{i}] reg={reg12[i]}  ref={ref['c12'][i]}  match={reg12[i]==ref['c12'][i]}")

# Reversal symmetry check for W11 coords at a generic point
import random
random.seed(1)
ptA = tuple(F(random.randint(1,9), random.randint(1,9)) for _ in range(5))
ptArev = ptA[::-1]
cA = coords_reg(ptA,11)
cArev = coords_reg(ptArev,11)
print("\nReversal check W11 at", ptA, "vs reversed:")
for i in (7,8):
    print(f"  [{i}] fwd={cA[i]}  rev={cArev[i]}  equal={cA[i]==cArev[i]}")
