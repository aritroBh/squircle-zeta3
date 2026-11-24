import mpmath as mp
mp.mp.dps = 60
s5=mp.sqrt(5); phi=(1+s5)/2; lf=mp.log(phi)
x1=-1/phi; x4=1/phi**2; kap1=1/phi; kap4=-1/phi**2
u1=lambda t: 1+kap1*t; u4=lambda t: 1+kap4*t
w1=lambda t: mp.log(u1(t)); w4=lambda t: mp.log(u4(t))
Li2=lambda x: mp.polylog(2,x)
def Li21(x,y):
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1); m=0
    while True:
        m+=1; xp*=x
        if m>=2: tot += xp*ps/m**2
        yp*=y; ps += yp/m
        if m>10 and abs(xp*ps/m**2) < mp.mpf('1e-70'): break
    return tot
Lam = Li21(x1,x1); swap = Li21(x4,-phi)
D3q = mp.quad(lambda t: mp.log(1-t)*w4(t)/u1(t), [0,1])
D4q = mp.quad(lambda t: mp.log(1-t)*w1(t)/u4(t), [0,1])
print("D3 quad:", mp.nstr(D3q,30))
print(" DEC1 form phi(lf*Li2(x4)+swap):", mp.nstr(phi*(lf*Li2(x4)+swap),30))
print(" my form   phi*swap-phi*lf*Li2(x4):", mp.nstr(phi*swap-phi*lf*Li2(x4),30))
print("D4 quad:", mp.nstr(D4q,30))
print(" DEC1 form phi^2(lf*Li2(x1)-Lam):", mp.nstr(phi**2*(lf*Li2(x1)-Lam),30))
# my reflection derivation for D4: D4 = int log(1-t) w1/u4 ; sub s=1-t
# 1+kap1(1-s) = (phi^2-s)/phi ; 1+kap4(1-s)= (phi+s)/phi^2 -> w1(1-s)= log(phi^2-s)-? 
# log(1+kap1(1-s)) = log(phi^2-s) - lf ; 1/u4(1-s) = phi/(phi+s)= 1/(1-(-1/phi)s)?? 
D1q = mp.quad(lambda t: mp.log(t)*w4(t)/u1(t), [0,1])
D2q = mp.quad(lambda t: mp.log(t)*w1(t)/u4(t), [0,1])
print("D1 quad:", mp.nstr(D1q,25), " Lam/x1:", mp.nstr(Lam/x1,25))
print("D2 quad:", mp.nstr(D2q,25), " swap/x4:", mp.nstr(swap/x4,25))
# candidate for D4 via reflection of D1: D4 = int log s * [w4(s)+?]/(phi(1+s/phi))? compute:
# 1+kap1(1-s) = (phi^2-s)/phi = phi*(1-s/phi^2)*... no: phi^2-s = phi^2(1-s/phi^2) -> (phi^2-s)/phi = phi*u4(s)
# log(1+kap1(1-s)) = log(phi^2-s)-lf = lf + log(1-s/phi^2) = lf + w4(s)
# 1+kap4(1-s) = (phi+s)/phi^2 = (1/phi)*(1+s/phi) -> u4(1-s) = (1/phi)*(1+s/phi) = u1(s)/phi
# so D4 = int log s (lf+w4(s)) * phi/u1(s) ds = phi*( lf*K1 + D1q ), K1=int log s/u1 ds
K1 = mp.quad(lambda t: mp.log(t)/u1(t), [0,1])
print("K1 =", mp.nstr(K1,25), " -Li2(x1)/kap1 =", mp.nstr(-Li2(x1)/kap1,25))
print("phi*(lf*K1+D1) =", mp.nstr(phi*(lf*K1+D1q),30), " vs D4 quad")
