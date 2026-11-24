import mpmath as mp
mp.mp.dps = 60
sq2 = mp.sqrt(2)
x1 = 1-sq2; x2 = 1-1/sq2
N = 300
def dsum(f):
    s = mp.mpf(0)
    for k1 in range(1,N):
        p1 = x1**k1
        for k2 in range(1,N):
            s += f(k1,k2)*p1*(x2**k2)
    return s
T = dsum(lambda k1,k2: 1/(k1*(k1+k2)))
Tp = dsum(lambda k1,k2: 1/(k2*(k1+k2)))
S1 = dsum(lambda k1,k2: (1/k1+1/k2)/(k1+k2))
print("T:", mp.nstr(T,30)); print("T':", mp.nstr(Tp,30)); print("T+T' vs S1:", mp.nstr(T+Tp-S1,3))
# So S1 = 2T. Compare with earlier: I111 = 2T - 2U fails, so maybe U part wrong:
U = dsum(lambda k1,k2: 1/(k1*(k1+k2+1)))
Up = dsum(lambda k1,k2: 1/(k2*(k1+k2+1)))
S2 = dsum(lambda k1,k2: (1/k1+1/k2)/(k1+k2+1))
print("U:", mp.nstr(U,30)); print("U':", mp.nstr(Up,30)); print("U+U' vs S2:", mp.nstr(U+Up-S2,3))
I111_sum = dsum(lambda k1,k2: 1/(k1*k2*(k1+k2+1)))
print("I111:", mp.nstr(I111_sum,30))
print("S1-S2:", mp.nstr(S1-S2,30))
print("2T-2U:", mp.nstr(2*T-2*U,30))
# test single term (1/k1+1/k2)/(k1+k2+1) vs my U decomposition pieces
k1,k2 = 2,3
lhs = (1/k1+1/k2)/(k1+k2+1)
print("single term lhs:", lhs)
