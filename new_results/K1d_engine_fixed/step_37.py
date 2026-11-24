import mpmath as mp
mp.mp.dps = 60
sq2 = mp.sqrt(2)
x1 = 1-sq2; x2 = 1-1/sq2
N = 400
def li11(a,b):
    s = mp.mpf(0)
    for n in range(2,N):
        an = a**n
        for m in range(1,n):
            s += an*(b**m)/(n*m)
    return s
l1 = -mp.log(1-x1); l2_ = -mp.log(1-x2)
T  = li11(x2, x1/x2)   # = Li_{1,1}(x2, x1/x2)
Tp = li11(x1, x2/x1)   # = Li_{1,1}(x1, x2/x1)
I111 = mp.quad(lambda t: mp.log(1+(1-x1)*t)*mp.log(1+(1-x2)*t), [0,1])  # kap = -x? NO
# careful: kappa_i = -x_i. kap1 = sqrt2-1 = -x1, kap2 = 1/sqrt2-1 = -x2
kap1, kap2 = -x1, -x2
I111 = mp.quad(lambda t: mp.log(1+kap1*t)*mp.log(1+kap2*t), [0,1])
formula = T*(1+1/x1) + Tp*(1+1/x2) - l1*(l2_-x2)/x2 - l2_*(l1-x1)/x1 - (l1-x1)/x1 - (l2_-x2)/x2
print("I111 direct :", mp.nstr(I111,40))
print("I111 formula:", mp.nstr(formula,40))
print("diff:", mp.nstr(I111-formula,3))
# stuffle partner check: T and T' independent? Li11(a,b)+Li11(b,a)=Li1(a)Li1(b)-Li2(ab)
a,b = x2, x1/x2
print("stuffle residual:", mp.nstr(li11(a,b)+li11(b,a) - (l1_:= -mp.log(1-a))*(-mp.log(1-b)) + mp.polylog(2,a*b),3))
