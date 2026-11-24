
from fractions import Fraction as F
import mzv_exact as M
import engine5, engine as E4

# ---- m0 part at a rational point: coords in DS basis ----
def gamprod_series_point(coeffs, Wmax):
    # coeffs: list of (ell_value (Fraction), exponent e=±1) for Gamma(1+ell*u)^e ; balanced
    # returns zmon series: list of dict sorted-tuple -> Fraction
    q=[{}]*(Wmax+1)
    for m in range(2,Wmax+1):
        s=sum(F(e)*ell**m for ell,e in coeffs)
        q[m]={(m,): F((-1)**m,m)*s}
    E=[{} for _ in range(Wmax+1)]; E[0]={():F(1)}
    for W in range(1,Wmax+1):
        acc={}
        for i in range(2,W+1):
            for z,c in q[i].items():
                for z2,c2 in E[W-i].items():
                    z3=tuple(sorted(z+z2))
                    acc[z3]=acc.get(z3,F(0))+F(i,W)*c*c2
        E[W]=acc
    return E

def m0_coords(a, Wmax):
    a1,a2,a3,a4,a5 = [F(x) for x in a]
    S4 = a1+a2+a3+a4; S4r = a5+a4+a3+a2
    s1 = gamprod_series_point([(a2+a4,1),(S4,-1),(a1+a3,1)], Wmax)   # * a2/(a2+a4)
    s2 = gamprod_series_point([(a2+a4,1),(S4r,-1),(a5+a3,1)], Wmax)  # * a4/(a2+a4)
    f1 = a2/(a2+a4); f2 = a4/(a2+a4)
    out={}
    for W in range(2,Wmax+1):
        T=engine5.table(W); d=T.nullity
        acc=[F(0)]*d
        for zmon,c in s1[W].items():
            c=c*f1
            wc=M.wordcombo(tuple((k,) for k in zmon))
            for w,cc in wc.items():
                rw=T.reduce_word(w)
                for b in range(d): acc[b]+=c*cc*rw[b]
        for zmon,c in s2[W].items():
            c=c*f2
            wc=M.wordcombo(tuple((k,) for k in zmon))
            for w,cc in wc.items():
                rw=T.reduce_word(w)
                for b in range(d): acc[b]+=c*cc*rw[b]
        out[W]=acc
    return out

def peval5(p, a):
    return sum(c*a[0]**e[0]*a[1]**e[1]*a[2]**e[2]*a[3]**e[3]*a[4]**e[4] for e,c in p.items())

forms6 = engine5.reduced_forms(6)[0]
# Control 1: a1=0 degeneration vs CHN1 P4 (a2,a3,a4,a5), W=2..6, several points
print("=== Control u1=0 (P5 regular + m0 == CHN1 P4) ===")
ok=True
for pt in [(1,1,1,1),(2,1,1,1),(1,2,3,2),(3,1,2,1)]:
    a5_=(0,)+pt
    for W in range(2,7):
        T=engine5.table(W); d=T.nullity
        mine=[peval5(forms6[W][b], a5_) if W>=3 else F(0) for b in range(d)]
        if W>=2:
            m0=m0_coords(a5_, W)[W]
            mine=[mine[b]+m0[b] for b in range(d)]
        ref=E4.reduced_coeffs(pt, W)[W]
        if mine!=ref:
            ok=False
            print("MISMATCH", pt, W, mine, ref)
print("u1=0 control:", "PASS" if ok else "FAIL")
