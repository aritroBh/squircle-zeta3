
import pickle, time
from fractions import Fraction as F
import mzv_exact as M
import engine as E4
import engine5, kimilib

forms10 = pickle.load(open('forms10.pkl','rb'))
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
def full_at(a, W):
    T=engine5.table(W); d=T.nullity
    reg=[kimilib.peval5(forms10[W][b], a) for b in range(d)]
    m0=kimilib.m0_coords(a, W)[W]
    return [reg[b]+m0[b] for b in range(d)]

t=time.time(); ok=True
for pt in [(1,1,1,1,1),(2,1,1,3,1),(1,2,3,2,1)]:
    a1,a2,a3,a4,a5 = pt
    for W in (9,10):
        if full_at((0,a2,a3,a4,a5), W) != E4.reduced_coeffs((a2,a3,a4,a5), W)[W]:
            ok=False; print("u1=0 FAIL", pt, W)
        if full_at((a1,a2,a3,a4,0), W) != E4.reduced_coeffs((a1,a2,a3,a4), W)[W]:
            ok=False; print("u5=0 FAIL", pt, W)
        ref3 = gamprod_coords([(a1,1),(a2,1),(a1+a2,-1),(a4,1),(a5,1),(a4+a5,-1)], W)
        if full_at((a1,a2,0,a4,a5), W) != ref3:
            ok=False; print("u3=0 FAIL", pt, W)
        ref2 = gamprod_coords([(a4,1),(a3+a5,1),(a3+a4+a5,-1)], W)
        if full_at((a1,0,a3,a4,a5), W) != ref2:
            ok=False; print("u2=0 FAIL", pt, W)
        ref4 = gamprod_coords([(a2,1),(a1+a3,1),(a1+a2+a3,-1)], W)
        if full_at((a1,a2,a3,0,a5), W) != ref4:
            ok=False; print("u4=0 FAIL", pt, W)
print("W=9,10 all degeneration controls:", "PASS" if ok else "FAIL", f"{time.time()-t:.0f}s")
revok=True
for W in (9,10):
    for b,p in enumerate(forms10[W]):
        pr = {(e[4],e[3],e[2],e[1],e[0]): c for e,c in p.items()}
        if pr!=p: revok=False; print("rev FAIL", W, b)
print("reversal symbolic W=9,10:", "PASS" if revok else "FAIL")
