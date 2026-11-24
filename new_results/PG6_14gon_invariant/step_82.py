import mpmath as mp, sys
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
mp.mp.dps=40
g=S.init()
Xv={i:E.fev(E.X[i],40) for i in E.IDX}
def Li21(x,y):
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1)
    for m in range(1,4000):
        xp*=x
        if m>=2: tot+=xp*ps/m**2
        yp*=y; ps+=yp/m
    return tot
def Gword(word):
    def rec(w,z):
        if not w: return mp.mpf(1)
        return mp.quad(lambda t: rec(w[1:],t)/(t-w[0]),[0,z],maxdegree=6)
    return rec(word,mp.mpf(1))
# verify G-word for more pairs
for (i,j) in [(2,5),(6,1),(5,6)]:
    lv=Li21(Xv[i],Xv[j]/Xv[i]); gw=Gword([mp.mpf(0),1/Xv[i],1/Xv[j]])
    print(f"Lam({i},{j}): diff {mp.nstr(abs(lv-gw),3)}")
# sigma(Lam) via S3(0, 1/x_i, 1/x_j) vs S_Li21
def ssub(A,B):
    R=dict(A)
    for k,c in B.items():
        R[k]=R.get(k,Fr(0))-c
        if R[k]==0: del R[k]
    return R
for (i,j) in [(1,2),(2,5)]:
    ri,rj = E.finv(E.X[i]), E.finv(E.X[j])
    s3 = S.S3(E.ZERO, ri, rj)
    sl21,_ = S.sym_Lam(i,j)
    d = ssub(s3, sl21)
    print(f"sigma Lam({i},{j}): S3(0,ri,rj) terms={len(s3)}, S_Li21 terms={len(sl21)}, diff terms={len(d)}")
    if 0<len(d)<20: print("   diff:", d)
