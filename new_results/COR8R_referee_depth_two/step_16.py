mp.mp.dps = 40
x, y = mp.mpf('0.3'), mp.mpf('0.4')
a = 1/y; b = 1/(x*y)

# numeric partial derivatives of G(a(x,y), b(x,y); 1) w.r.t. log x and log y
h = mp.mpf('1e-12')
def G_of(xx, yy):
    return G_int([1/yy, 1/(xx*yy)], 1)
dlogx = (G_of(x*mp.exp(h), y) - G_of(x*mp.exp(-h), y))/(2*h)
dlogy = (G_of(x, y*mp.exp(h)) - G_of(x, y*mp.exp(-h)))/(2*h)

def Li01(xx, yy):  # Li_{0,1}(x,y) = sum_{k1<k2} x^k1 y^k2 / k2
    s = mp.mpf(0)
    for k2 in range(2, 300):
        inner = sum(xx**k1 for k1 in range(1, k2))
        s += yy**k2/k2 * inner
    return s
def Li10(xx, yy):  # Li_{1,0}(x,y) = sum_{k1<k2} x^k1 y^k2 / k1
    s = mp.mpf(0)
    for k2 in range(2, 300):
        inner = sum(xx**k1/k1 for k1 in range(1, k2))
        s += yy**k2 * inner
    return s

print("dG/dlogx numeric:", mp.nstr(dlogx, 25))
print("Li_{0,1}(x,y)   :", mp.nstr(Li01(x,y), 25))
print("dG/dlogy numeric:", mp.nstr(dlogy, 25))
print("Li_{1,0}(x,y)   :", mp.nstr(Li10(x,y), 25))
