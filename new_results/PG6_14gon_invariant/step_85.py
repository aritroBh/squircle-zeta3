import mpmath as mp, sys
sys.path.insert(0,'.')
import engine14 as E
mp.mp.dps=40
Xv={i:E.fev(E.X[i],40) for i in E.IDX}
# prefix f(x,y) = sum_{0<k<m} x^m y^k/(m^2 k) = Li_{1,2}(y,x)
def prefix(x,y):
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
x=Xv[1]; y=Xv[2]/Xv[1]; xy=x*y
pv = prefix(x,y)   # = Li_{1,2}(y,x)
print("prefix f(x,y)=Li_{1,2}(y,x) =", mp.nstr(pv,14))
# candidate G-words for Li_{1,2}(y,x)
for nm,wd in {'I(1/(xy),0,1/x)':[1/xy,mp.mpf(0),1/x],'I(1/(xy),0,1/y)':[1/xy,mp.mpf(0),1/y],'I(0,1/(xy),1/x)':[mp.mpf(0),1/xy,1/x],'I(0,1/x,1/(xy))':[mp.mpf(0),1/x,1/xy],'I(1/y,0,1/(xy))':[1/y,mp.mpf(0),1/xy]}.items():
    try:
        v=Gword(wd); print(f"  {nm}: {mp.nstr(v,14)} diff {mp.nstr(abs(pv-v),2)}")
    except Exception as ex: print(f"  {nm}: FAIL {str(ex)[:40]}")
