import mpmath as mp, sys
sys.path.insert(0,'.')
import engine14 as E
mp.mp.dps=40
Xv={i:E.fev(E.X[i],40) for i in E.IDX}
def Li21(x,y):
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1)
    for m in range(1,4000):
        xp*=x
        if m>=2: tot+=xp*ps/m**2
        yp*=y; ps+=yp/m
    return tot
# fast evaluator for G(0,a,b;1) and G(a,b,c;1) using series for the innermost + quad
def Gword(word):
    # evaluate iterated integral via recursive quad but with low maxdegree for speed
    def rec(w, z):
        if not w: return mp.mpf(1)
        a=w[0]
        return mp.quad(lambda t: rec(w[1:], t)/(t-a), [0,z], maxdegree=6)
    return rec(word, mp.mpf(1))
x1,x2=Xv[1],Xv[2]; y=x2/x1
lv=Li21(x1,x2/x1)
print("Lam(1,2)=",mp.nstr(lv,14))
cands = {
 'G(0,1/x2,x1/x2)':[mp.mpf(0),1/x2,x1/x2],
 '-G(0,1/x2,x1/x2)':'NEG0',
 'G(0,1/x2,1/x1)':[mp.mpf(0),1/x2,1/x1],
 'G(1/x1,0,1/x2)':[1/x1,mp.mpf(0),1/x2],
 'G(0,1/x1,1/x2)':[mp.mpf(0),1/x1,1/x2],
}
for nm,wd in cands.items():
    if wd=='NEG0': continue
    try:
        v=Gword(wd)
        print(f"  {nm}: {mp.nstr(v,14)} diff {mp.nstr(abs(lv-v),2)} absratio {mp.nstr(abs(lv/v) if v else 0,4)}")
    except Exception as ex:
        print(f"  {nm}: FAIL {str(ex)[:50]}")
