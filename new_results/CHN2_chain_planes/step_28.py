
import mpmath as mp
mp.mp.dps = 50
u1,u2,u3,u4,u5 = mp.mpf('0.5'), mp.mpf(1)/3, mp.mpf('0.25'), mp.mpf('0.2'), mp.mpf(1)/6
def I2(T,a1,a2):
    return T**a2 * mp.hyper([a2,-a1],[1+a2], T)
for t3 in [mp.mpf('0.5'), mp.mpf('0.01'), mp.mpf('1e-10'), mp.mpf('1e-40')]:
    try:
        v = I2(1-t3, u1, u2)
        print(t3, "->", mp.nstr(v, 20))
    except Exception as e:
        print(t3, "ERROR", e)
# hyper at z=1 exactly?
print("hyper at 1:", mp.nstr(mp.hyper([u2,-u1],[1+u2], 1), 20))
print("gamma form:", mp.nstr(mp.gamma(1+u2)*mp.gamma(1+u1)/mp.gamma(1+u1+u2), 20))
