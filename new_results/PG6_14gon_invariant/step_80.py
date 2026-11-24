import mpmath as mp, sys
sys.path.insert(0,'.')
import engine14 as E
mp.mp.dps=60
Xv={i:E.fev(E.X[i],60) for i in E.IDX}
def Li21(x,y):
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1)
    for m in range(1,6000):
        xp*=x
        if m>=2: tot+=xp*ps/m**2
        yp*=y; ps+=yp/m
    return tot
# generic G(word;1) evaluator, safe for 0 letters (inner integral vanishes at 0)
def Geval(word):
    if not word: return mp.mpf(1)
    a=word[0]
    return mp.quad(lambda t: Geval(word[1:]) if False else _G(word[1:],t)/(t-a), [0,1])
def _G(word,t):
    if not word: return mp.mpf(1)
    a=word[0]
    # int_0^t ds/(s-a) * _G(word[1:], s)
    return mp.quad(lambda s: _G(word[1:],s)/(s-a), [0,t])
# candidate G-word forms for Lam(1,2)=Li21(x1, x2/x1);  letters to try
x1,x2=Xv[1],Xv[2]; y=x2/x1
lv=Li21(x1,x2/x1)
cands = {
 'I(1/(xy),0,1/x)':[1/(x1*y),0,1/x1],
 'I(1/x,0,1/(xy))':[1/x1,0,1/(x1*y)],
 'I(1/y,0,1/(xy))':[1/y,0,1/(x1*y)],
 'I(1/(xy),0,1/y)':[1/(x1*y),0,1/y],
 'I(1/x,1/(xy),0)':[1/x1,1/(x1*y),0],
 'I(0,1/(xy),1/x)':[0,1/(x1*y),1/x1],
 'I(1/x,0,0)':[1/x1,0,0],
}
print("Lam(1,2) =", mp.nstr(lv,16))
for nm,wd in cands.items():
    try:
        v=Geval(wd)
        print(f"  {nm}: {mp.nstr(v,16)}  diff {mp.nstr(abs(lv-v),2)}  ratio {mp.nstr(lv/v if v else 0,4)}")
    except Exception as ex:
        print(f"  {nm}: FAIL {str(ex)[:40]}")
