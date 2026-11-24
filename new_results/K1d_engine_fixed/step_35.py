import mpmath as mp
mp.mp.dps = 60
sq2 = mp.sqrt(2)
kap1 = sq2-1; kap2 = 1/sq2-1
x1 = -kap1; x2 = -kap2
N = 300
def dsum(f):
    s = mp.mpf(0)
    for k1 in range(1,N):
        p1 = x1**k1
        for k2 in range(1,N):
            s += f(k1,k2)*p1*(x2**k2)
    return s
I111_sum = dsum(lambda k1,k2: 1/(k1*k2*(k1+k2+1)))
I111 = mp.quad(lambda t: mp.log(1+kap1*t)*mp.log(1+kap2*t), [0,1])
print("I111 quad:", mp.nstr(I111,40))
print("I111 sum :", mp.nstr(I111_sum,40), " diff:", mp.nstr(I111-I111_sum,3))
T = dsum(lambda k1,k2: 1/(k1*(k1+k2)))
U = dsum(lambda k1,k2: 1/(k1*(k1+k2+1)))
print("2T-2U =", mp.nstr(2*T-2*U,40), " diff vs I111:", mp.nstr(2*T-2*U-I111,3))
U1 = dsum(lambda k1,k2: 1/(k1*(k2+1)))
U2 = dsum(lambda k1,k2: 1/((k1+k2+1)*(k2+1)))
print("U vs U1-U2:", mp.nstr(U-(U1-U2),3))
l1 = -mp.log(1-x1); l2_ = -mp.log(1-x2)
U1f = l1*(l2_-x2)/x2
print("U1:", mp.nstr(U1,30), " formula:", mp.nstr(U1f,30), " diff:", mp.nstr(U1-U1f,3))
# U2 check: U2 = (1/x2)[Li11(x1,x2/x1) - x2*(l1-x1)/x1]  -- verify Li11 numerically by truncation
def li11(a,b):
    s = mp.mpf(0)
    for n in range(2,N):
        an = a**n
        for m in range(1,n):
            s += an*(b**m)/(n*m)
    return s
Bv = li11(x1, x2/x1)
m1part = x2*(l1-x1)/x1
U2f = (Bv - m1part)/x2
print("U2:", mp.nstr(U2,30), " formula:", mp.nstr(U2f,30), " diff:", mp.nstr(U2-U2f,3))
