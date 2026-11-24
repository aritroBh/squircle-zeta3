import mpmath as mp, sys
sys.path.insert(0,'.')
import engine14 as E
mp.mp.dps=60
Xv={i: E.fev(E.X[i],60) for i in E.IDX}
def Geval(word, z=1):
    if not word: return mp.mpf(1)
    a=word[0]; rest=word[1:]
    return mp.quad(lambda t: Geval(rest, t)/(t-a), [0,z])
w={i:(lambda tt,i=i: mp.log(1-Xv[i]*tt)) for i in E.IDX}
u={i:(lambda tt,i=i: 1-Xv[i]*tt) for i in E.IDX}
r={i: 1/Xv[i] for i in E.IDX}
# T2(1,2) = -(1/x1) G(1/x1,1/x2)
v=mp.quad(lambda tt: w[2](tt)/u[1](tt),[0,1]); g=-r[1]*Geval([r[1],r[2]])
print("T2(1,2):", mp.nstr(v,14), mp.nstr(g,14), "diff", mp.nstr(abs(v-g),3))
# T3a(1,2,5) = -(1/x1)[G(r1,r2,r5)+G(r1,r5,r2)]
v=mp.quad(lambda tt: w[2](tt)*w[5](tt)/u[1](tt),[0,1]); g=-r[1]*(Geval([r[1],r[2],r[5]])+Geval([r[1],r[5],r[2]]))
print("T3a(1,2,5):", mp.nstr(v,14), mp.nstr(g,14), "diff", mp.nstr(abs(v-g),3))
# T3p(1,2) = -(2/x1) G(r1,r2,r2)
v=mp.quad(lambda tt: w[2](tt)**2/u[1](tt),[0,1]); g=-2*r[1]*Geval([r[1],r[2],r[2]])
print("T3p(1,2):", mp.nstr(v,14), mp.nstr(g,14), "diff", mp.nstr(abs(v-g),3))
# Lam(1,2) = Li_{2,1}(x1, x2/x1): G-word? Li_{2,1}(x,y)=? try I(0; 1/(xy),1/y;1)? 
def Li21(x,y):
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1)
    for m in range(1,6000):
        xp*=x
        if m>=2: tot+=xp*ps/m**2
        yp*=y; ps+=yp/m
    return tot
lv = Li21(Xv[1], Xv[2]/Xv[1])
# candidate: Li_{2,1}(x,y) = -G(1/(xy), 1/x ... let me test Li_{2,1}(x,y) = -G(1/y, 1/(xy);1)? 
for cand,dsc in [([1/(Xv[2]/Xv[1]), 1/Xv[1]],'G(1/y,1/x)'),([1/Xv[1],1/(Xv[2]/Xv[1])],'G(1/x,1/y)')]:
    print("Lam(1,2)=",mp.nstr(lv,14)," cand",dsc,"=",mp.nstr(-Geval(cand),14), mp.nstr(Geval(cand),14))
