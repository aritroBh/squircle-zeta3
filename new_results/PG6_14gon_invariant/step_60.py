import mpmath as mp, sys
sys.path.insert(0,'.')
import engine14 as E
mp.mp.dps=60
Xv={i: E.fev(E.X[i],60) for i in E.IDX}
def Geval(word, z=1):
    if not word: return mp.mpf(1)
    a=word[0]; rest=word[1:]
    return mp.quad(lambda t: Geval(rest, t)/(t-a), [0,z])
def Li21(x,y):
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1)
    for m in range(1,6000):
        xp*=x
        if m>=2: tot+=xp*ps/m**2
        yp*=y; ps+=yp/m
    return tot
def Aser(x):
    tot=mp.mpf(0); xp=mp.mpf(1); H=mp.mpf(0)
    for n in range(1,8000): xp*=x; H+=mp.mpf(1)/n; tot+=H*xp/n**2
    return tot
r={i: 1/Xv[i] for i in E.IDX}
# Lam(1,2) = Li21(x1, x2/x1). Try I(0; r2, 0, r1; 1) [Li_{1,2}(y,x) form]
lv = Li21(Xv[1], Xv[2]/Xv[1])
for word,dsc in [([r[2],0,r[1]],'I(r2,0,r1)'),([0,r[2],r[1]],'I(0,r2,r1)'),([r[1],0,r[2]],'I(r1,0,r2)')]:
    print("Lam(1,2)=",mp.nstr(lv,14), dsc,"=",mp.nstr(Geval(word),14)," diff",mp.nstr(abs(lv-Geval(word)),3))
# Li3(x1): I(0; r1,0,0;1) = -Li3? Goncharov Li_n(z) = -I(0;1/z,0..0;1)? test
lv3 = mp.polylog(3,Xv[1])
for word,dsc in [([r[1],0,0],'I(r1,0,0)'),([0,0,r[1]],'I(0,0,r1)')]:
    print("Li3(x1)=",mp.nstr(lv3,14), dsc,"=",mp.nstr(Geval(word),14), " -that=",mp.nstr(-Geval(word),14))
# Ax(1) = A(x1) = Li3(x1) + (1/2) int_0^x1 log^2(1-t)/t dt
av = Aser(Xv[1])
# int_0^x log^2(1-t)/t dt = ? G-word: = 2*I(0; 1/x1 ... let me just compute Aser symbol later via Li3 + repeated repair
print("Ax(1)=", mp.nstr(av,14), " Li3 part=", mp.nstr(lv3,14), " diff=", mp.nstr(av-lv3,14))
open('progress.log','a').write('PG6 s_B7: G-word forms for T2/T3a/T3p verified; Lam/Li3 forms being pinned\n')
