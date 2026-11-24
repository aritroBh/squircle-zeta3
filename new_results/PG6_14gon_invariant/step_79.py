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
# efficient G(a,0,c;1) via nested quad with analytic inner level
def G_a0c(a,c):
    # G(a,0,c;1) = int_0^1 dt/(t-a) * G(0,c;t); G(0,c;t)= int_0^t ds/s * log((s-c)/(-c))... careful
    def G0c(t):
        # int_0^t (1/s) * [int_0^s du/(u-c)] ds ; inner = log((s-c)/(-c))
        return mp.quad(lambda s: mp.log((s-c)/(-c))/s, [0,t])
    return mp.quad(lambda t: G0c(t)/(t-a), [0,1])
for (i,j) in [(1,2),(2,5)]:
    lv = Li21(Xv[i], Xv[j]/Xv[i])
    gw = G_a0c(1/Xv[j], 1/Xv[i])   # I(0; 1/x_j, 0, 1/x_i; 1) = G(1/x_j, 0, 1/x_i;1)
    print(f"Lam({i},{j}) = {mp.nstr(lv,16)}  G-word = {mp.nstr(gw,16)}  diff {mp.nstr(abs(lv-gw),3)}")
