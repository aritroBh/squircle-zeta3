import mpmath as mp
mp.mp.dps = 25
c0 = mp.mpf('0.31'); B0 = 2/(1+c0)
# numeric: Li_{1,1,1}(x,y,z) = sum_{k1<k2<k3} x^k1 y^k2 z^k3/(k1 k2 k3)
def Li111(x,y,z,K=120):
    s = mp.mpf(0)
    for k3 in range(3,K):
        for k2 in range(2,k3):
            for k1 in range(1,k2):
                s += x**k1*y**k2*z**k3/(k1*k2*k3)
    return s
val_Li = Li111(1/B0, 2, mp.mpf('0.5'))
val_G20B = mp.quad(lambda t1: mp.quad(lambda t2: mp.quad(lambda t3: 1/(t3-B0), [mp.mpf(0),t2])/(t2), [mp.mpf(0),t1])/(t1-2), [mp.mpf(0),mp.mpf(1)])
print("G(2,0,B;1)        =", mp.nstr(val_G20B,18))
print("-Li_{1,1,1}(1/B,2,1/2) =", mp.nstr(-val_Li,18))
# and G(2,1,B;1)?
val_G21B = mp.quad(lambda t1: mp.quad(lambda t2: mp.quad(lambda t3: 1/(t3-B0), [mp.mpf(0),t2])/(t2-1), [mp.mpf(0),t1])/(t1-2), [mp.mpf(0),mp.mpf(1)])
print("G(2,1,B;1)        =", mp.nstr(val_G21B,18))
print("Li_{1,1,1}(1/B,2,1/2)  =", mp.nstr(val_Li,18))
