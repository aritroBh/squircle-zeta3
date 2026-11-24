
from fractions import Fraction as F
import mzv_exact as M
import engine as E4
import engine5, kimilib

def gamprod_coords(coeffs, W):
    s = kimilib.gamprod_series_point(coeffs, W)
    T=engine5.table(W); d=T.nullity
    acc=[F(0)]*d
    for zmon,c in s[W].items():
        wc=M.wordcombo(tuple((k,) for k in zmon))
        for w,cc in wc.items():
            rw=T.reduce_word(w)
            for b in range(d): acc[b]+=c*cc*rw[b]
    return acc

forms6, eng6 = engine5.reduced_forms(6)
def full_at(a, W):
    T=engine5.table(W); d=T.nullity
    reg=[kimilib.peval5(forms6[W][b], a) if W>=3 else F(0) for b in range(d)]
    m0=kimilib.m0_coords(a, W)[W]
    return [reg[b]+m0[b] for b in range(d)]

ok = {'u1=0':True,'u5=0':True,'u3=0':True,'u2=0':True,'u4=0':True}
pts5 = [(1,1,1,1,1),(2,1,1,3,1),(1,2,3,2,1),(3,1,2,1,2),(2,3,1,2,3)]
for pt in pts5:
    a1,a2,a3,a4,a5 = pt
    for W in range(2,7):
        # u1=0
        if full_at((0,a2,a3,a4,a5), W) != E4.reduced_coeffs((a2,a3,a4,a5), W)[W]:
            ok['u1=0']=False; print("u1=0 FAIL", pt, W)
        # u5=0
        if full_at((a1,a2,a3,a4,0), W) != E4.reduced_coeffs((a1,a2,a3,a4), W)[W]:
            ok['u5=0']=False; print("u5=0 FAIL", pt, W)
        # u3=0
        ref = gamprod_coords([(a1,1),(a2,1),(a1+a2,-1),(a4,1),(a5,1),(a4+a5,-1)], W)
        if full_at((a1,a2,0,a4,a5), W) != ref:
            ok['u3=0']=False; print("u3=0 FAIL", pt, W, full_at((a1,a2,0,a4,a5),W), ref)
        # u2=0
        ref = gamprod_coords([(a4,1),(a3+a5,1),(a3+a4+a5,-1)], W)
        if full_at((a1,0,a3,a4,a5), W) != ref:
            ok['u2=0']=False; print("u2=0 FAIL", pt, W, full_at((a1,0,a3,a4,a5),W), ref)
        # u4=0
        ref = gamprod_coords([(a2,1),(a1+a3,1),(a1+a2+a3,-1)], W)
        if full_at((a1,a2,a3,0,a5), W) != ref:
            ok['u4=0']=False; print("u4=0 FAIL", pt, W, full_at((a1,a2,a3,0,a5),W), ref)
print(ok)
# reversal invariance of symbolic forms W=3..6
revok=True
for W in range(3,7):
    for b,p in enumerate(forms6[W]):
        pr = {(e[4],e[3],e[2],e[1],e[0]): c for e,c in p.items()}
        if pr!=p: revok=False; print("rev FAIL", W, b)
print("reversal symbolic:", "PASS" if revok else "FAIL")
