import mpmath as mp
mp.mp.dps = 60

def Li21_series(x, y):
    # sum_{0<m<n} x^m y^n / (m^2 n), convergent for |x|<=1,|y|<=1 (not both on boundary badly)
    x = mp.mpf(x); y = mp.mpf(y)
    s = mp.mpf(0); N = 400
    yn = y
    for n in range(2, N):
        yn *= y  # y^n
        t = mp.mpf(0); xm = x
        for m in range(1, n):
            t += xm/m**2
            xm *= x
        s += t*yn/n
    return s

def Li21_int(x, y):
    # Li_{2,1}(x,y) = int_0^1 y Li_2(x y t)/(1 - y t) dt   (derived in analysis)
    x = mp.mpf(x); y = mp.mpf(y)
    f = lambda t: y*mp.polylog(2, x*y*t)/(1-y*t)
    return mp.quad(f, [0, 1])

def Li21_int_s(x,y):
    # int_0^y Li_2(x s)/(1-s) ds
    x = mp.mpf(x); y = mp.mpf(y)
    return mp.quad(lambda s: mp.polylog(2, x*s)/(1-s), [0, y])

for (x,y) in [(0.3,0.4), (-0.5,0.7), (0.2,-0.6), (-0.8,-0.3)]:
    a = Li21_series(x,y); b = Li21_int(x,y); c = Li21_int_s(x,y)
    print((x,y), "series-int:", mp.nstr(abs(a-b),3), " series-ints:", mp.nstr(abs(a-c),3))

# also check G(1/y, 0, 1/(xy); 1) explicitly with triple integral at one point
x, y = 0.3, 0.4
a1, a2, a3 = 1/y, 0.0, 1/(x*y)
G = mp.quad(lambda t1: 1/(t1-a1)*mp.quad(lambda t2: 1/(t2-a2)*mp.quad(lambda t3: 1/(t3-a3), [0,t2]), [0,t1]), [0,1])
print("G(1/y,0,1/(xy);1) vs series:", mp.nstr(abs(G - Li21_series(x,y)),3))
