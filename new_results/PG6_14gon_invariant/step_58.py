import mpmath as mp, sys
sys.path.insert(0,'.')
import engine14 as E
mp.mp.dps=60
Xv={i: E.fev(E.X[i],60) for i in E.IDX}
def Geval(word, z=1, pts=None):
    if not word: return mp.mpf(1)
    a=word[0]; rest=word[1:]
    if a==0:
        # dt/t singular at 0; Geval(rest,0)=0, integrand ~ finite. split near 0.
        return mp.quad(lambda t: Geval(rest,t)/(t-a) if t>1e-30 else mp.mpf(0), [mp.mpf('1e-30'),z]) + \
               mp.quad(lambda t: Geval(rest,t)/t, [0,mp.mpf('1e-30')])
    return mp.quad(lambda t: Geval(rest, t)/(t-a), [0,z])
def Li21(x,y):
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1)
    for m in range(1,6000):
        xp*=x
        if m>=2: tot+=xp*ps/m**2
        yp*=y; ps+=yp/m
    return tot
r={i: 1/Xv[i] for i in E.IDX}
# Li_{2,1}(x1, x2/x1): G(0, 1/x1, 1/x2 ; 1)
lv = Li21(Xv[1], Xv[2]/Xv[1])
g = Geval([0, r[1], r[2]])
print("Lam(1,2)=Li21(x1,x2/x1):", mp.nstr(lv,15), " G(0,r1,r2):", mp.nstr(g,15), "diff", mp.nstr(abs(lv-g),3))
# also with sign
# Li_3(z): = -G(0,0,1/z;1)? test z=x1
z=Xv[1]
li3 = mp.polylog(3,z)
g3a = -Geval([0,0,1/z]); g3b = Geval([0,0,1/z])
print("Li3(x1):", mp.nstr(li3,15), " -G(0,0,1/z):", mp.nstr(g3a,15), " +G:", mp.nstr(g3b,15))
# A(x) = Li3(x) + (1/2) int_0^x log^2(1-t)/t dt ; the latter = ? G-word: int_0^x log^2(1-t)/t dt = 2*G-word?
# log(1-t)=G(1;t). int_0^x [G(1;t)]^2/t dt = 2 G(0,1,1; x)? with endpoint x. Let me test at x=x1
x=Xv[1]
Aser=lambda xx: sum((mp.mpf(1)) for _ in [0]) # placeholder
def Aser(x):
    tot=mp.mpf(0); xp=mp.mpf(1); H=mp.mpf(0)
    for n in range(1,8000): xp*=x; H+=mp.mpf(1)/n; tot+=H*xp/n**2
    return tot
Av = Aser(x)
part2 = mp.quad(lambda t: mp.log(1-t)**2/t, [0,x])
print("A(x1):", mp.nstr(Av,15), " Li3 + 0.5*int:", mp.nstr(mp.polylog(3,x)+0.5*part2,15))
# int_0^x log^2(1-t)/t dt as G-word: = 2 G(0,1,1; x)? test
g2 = 2*Geval([0,1,1], x)
print("int log^2(1-t)/t:", mp.nstr(part2,15), " 2G(0,1,1;x):", mp.nstr(g2,15), "diff", mp.nstr(abs(part2-g2),3))
