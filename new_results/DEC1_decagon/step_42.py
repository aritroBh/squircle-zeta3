import os, mpmath as mp
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 100
s5=mp.sqrt(5); phi=(1+s5)/2; lf=mp.log(phi); z2=mp.pi**2/6
x1=-1/phi; x4=1/phi**2
kap1=1/phi; kap4=-1/phi**2
u1=lambda t: 1+kap1*t; u4=lambda t: 1+kap4*t
w1=lambda t: mp.log(u1(t)); w4=lambda t: mp.log(u4(t))
Li2=lambda x: mp.polylog(2,x); Li3=lambda x: mp.polylog(3,x)
def Li21(x,y):
    # sum_{m>n>=1} x^m y^n/(m^2 n): iterate m, keep partial sum of y^n/n
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1)
    m=0
    while True:
        m+=1
        xp*=x
        if m>=2:
            tot += xp*ps/m**2
        yp*=y
        ps += yp/m
        if m>10 and abs(xp*ps/m**2) < mp.mpf('1e-110'): break
    return tot,m
v,m = Li21(x1,x1)
print("Li21(x1,x1) series =", mp.nstr(v,30), "m=",m)
print("Li3(x1^2)=Li3(x4)   =", mp.nstr(Li3(x4),30))
print("ratio Li21/Li3:", mp.nstr(v/Li3(x4),15))
# stuffle & shuffle pieces
Li12_xx = None
def Li12(x,y):
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1); m=0
    while True:
        m+=1
        xp*=x
        if m>=2: tot += xp*ps/m
        yp*=y
        ps += yp/m**2
        if m>10 and abs(xp*ps/m) < mp.mpf('1e-110'): break
    return tot
v12 = Li12(x1,x1)
print("Li12(x1,x1) =", mp.nstr(v12,30))
L1x = -mp.log(1-x1); L2x = Li2(x1)
print("stuffle: L1*L2 - (Li21+Li12+Li3(x^2)) =", mp.nstr(L1x*L2x-(v+v12+Li3(x4)),8))
print("shuffle: L1*L2 - (Li12+2Li21) =", mp.nstr(L1x*L2x-(v12+2*v),8))
# D1 via series formula: Li21(x1, x4/x1)/x1
D1_series = v/x1
D1q = mp.quad(lambda t: mp.log(t)*w4(t)/u1(t),[0,1])
print("D1 quad:", mp.nstr(D1q,25), " D1 series:", mp.nstr(D1_series,25), " diff:", mp.nstr(abs(D1q-D1_series),6))
# reflection identity for D4 by quadrature both sides
D4q = mp.quad(lambda t: mp.log(1-t)*w1(t)/u4(t),[0,1])
D4ref = phi**2*mp.quad(lambda s: mp.log(s)*(mp.log(phi**2-s)-lf)/(phi+s),[0,1])
print("D4 quad:", mp.nstr(D4q,25), " D4 reflected quad:", mp.nstr(D4ref,25), " diff:", mp.nstr(abs(D4q-D4ref),6))
# D4 in terms of D1q:
K = mp.quad(lambda s: mp.log(s)/u1(s),[0,1])
print("K = -Li2(x1)/x1:", mp.nstr(K,25), mp.nstr(-Li2(x1)/x1,25))
D4c = phi*(lf*K + D1q)
print("D4 = phi*(lf*K+D1):", mp.nstr(D4c,25), " diff vs quad:", mp.nstr(abs(D4c-D4q),6))