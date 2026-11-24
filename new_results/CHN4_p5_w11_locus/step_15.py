import sys, pickle
sys.path.insert(0,'.')
import point_engine as PE, mzv_exact as M, kimilib
from fractions import Fraction as F

def coords_reg(a, W):
    vw = PE.volreg_words_point(a, W); vwr = PE.volreg_words_point(a[::-1], W)
    T = PE.table(W); d=T.nullity; acc=[F(0)]*d
    for src in (vw[W], vwr[W]):
        for w,pp in src.items():
            rw=T.reduce_word(w)
            for b in range(d):
                if rw[b]: acc[b]+=pp*rw[b]
    return acc

pt=(1,2,3,4,5)
ref = pickle.load(open('test_pt_w1112.pkl','rb'))
reg12 = coords_reg(pt,12)
m0_12 = kimilib.m0_coords(pt,12)[12]
full12 = [reg12[b]+m0_12[b] for b in range(12)]
# CLASS = idx9 + (14/9)*idx10
regCLASS = reg12[9] + F(14,9)*reg12[10]
fullCLASS = full12[9] + F(14,9)*full12[10]
refCLASS = ref['c12'][9] + F(14,9)*ref['c12'][10]
print("m0 contribution to idx9:", m0_12[9])
print("m0 contribution to idx10:", m0_12[10])
print("m0 contribution to CLASS:", m0_12[9]+F(14,9)*m0_12[10])
print()
print("regCLASS  =", regCLASS)
print("fullCLASS =", fullCLASS)
print("refCLASS  =", refCLASS)
print("reg==full CLASS:", regCLASS==fullCLASS, " full==ref CLASS:", fullCLASS==refCLASS)
print()
print("idx11 reg/full/ref:", reg12[11], full12[11], ref['c12'][11], "m0 idx11:", m0_12[11])
print("idx7 reg/full/ref:", reg12[7], full12[7], ref['c12'][7], "m0 idx7:", m0_12[7])
print("idx8 reg/full/ref:", reg12[8], full12[8], ref['c12'][8], "m0 idx8:", m0_12[8])
