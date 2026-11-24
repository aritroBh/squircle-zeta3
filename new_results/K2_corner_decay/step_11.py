
import mpmath as mp
mp.mp.dps = 60

def ts_sum(f, h, J):
    tot = mp.mpf(0)
    for jj in range(-J, J+1):
        t = jj*h; st = mp.sinh(t)
        g = mp.tanh(mp.pi/2*st)
        w = (mp.pi/2)*mp.cosh(t)/(mp.cosh(mp.pi/2*st)**2)
        tot += w*f(g)
    return tot*h

# test 1: int_0^1 1 dx = 1
for h,J in [(mp.mpf('0.16'),34),(mp.mpf('0.08'),68)]:
    print("int 1:", mp.nstr(ts_sum(lambda g: 1, h, J)-1, 5))
# test 2: int_0^1 sqrt(1+g^2)/(1+g) dx  (smooth)
f2 = lambda g: mp.sqrt(1+g*g)/(1+g)
v2 = mp.quad(f2, [0,1])
for h,J in [(mp.mpf('0.16'),34),(mp.mpf('0.08'),68),(mp.mpf('0.04'),136)]:
    print("int f2 err:", mp.nstr(ts_sum(f2, h, J)-v2, 5))
# test 3: int_0^1 g^6 dx = 1/7 (endpoint vanishing)
for h,J in [(mp.mpf('0.16'),34),(mp.mpf('0.08'),68)]:
    print("int g^6 err:", mp.nstr(ts_sum(lambda g: g**6, h, J)-mp.mpf(1)/7, 5))
