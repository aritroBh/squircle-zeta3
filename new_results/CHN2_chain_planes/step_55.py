
from fractions import Fraction as F
import mzv_exact as M
import engine as E4
import engine5, kimilib

def m0_words(a, Wmax):
    a1,a2,a3,a4,a5 = [F(x) for x in a]
    S4 = a1+a2+a3+a4; S4r = a5+a4+a3+a2
    s1 = kimilib.gamprod_series_point([(a2+a4,1),(S4,-1),(a1+a3,1)], Wmax)
    s2 = kimilib.gamprod_series_point([(a2+a4,1),(S4r,-1),(a5+a3,1)], Wmax)
    f1 = a2/(a2+a4); f2 = a4/(a2+a4)
    out={}
    for W in range(2,Wmax+1):
        acc={}
        for s_,f_ in ((s1,f1),(s2,f2)):
            for zmon,c in s_[W].items():
                c=c*f_
                wc=M.wordcombo(tuple((k,) for k in zmon))
                for w,cc in wc.items():
                    acc[w]=acc.get(w,F(0))+c*cc
        out[W]={w:c for w,c in acc.items() if c}
    return out

eng6 = engine5.build_engine(6)
vw = engine5.volreg_words(eng6, 6)
pt=(1,1,1,1); a5_=(0,)+pt
W=5
tw = E4.coeffs(pt, W)[W]
m0w = m0_words(a5_, W)[W]
treg = dict(tw)
for w,c in m0w.items():
    treg[w]=treg.get(w,F(0))-c
    if treg[w]==0: del treg[w]
mine = {w: kimilib.peval5(p, a5_) for w,p in vw[W].items()}
allw = sorted(set(treg)|set(mine))
for w in allw:
    t=treg.get(w,F(0)); m=mine.get(w,F(0))
    flag = "" if t==m else "   <<<< DIFF"
    print(f"  {w}: truth={t} mine={m}{flag}")
