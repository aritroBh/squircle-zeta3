from cor8r_num import *
setdps(40)
x, y = mp.mpf('0.3'), mp.mpf('0.4')
h = mp.mpf('1e-12')
def G_of(xx, yy):
    return G_int([1/yy, 1/(xx*yy)], 1)
dlogx = (G_of(x*mp.exp(h), y) - G_of(x*mp.exp(-h), y))/(2*h)
dlogy = (G_of(x, y*mp.exp(h)) - G_of(x, y*mp.exp(-h)))/(2*h)

def Li01(xx, yy):
    s = mp.mpf(0)
    for k2 in range(2, 400):
        s += yy**k2/k2 * sum(xx**k1 for k1 in range(1, k2))
    return s
def Li10(xx, yy):
    s = mp.mpf(0)
    for k2 in range(2, 400):
        s += yy**k2 * sum(xx**k1/k1 for k1 in range(1, k2))
    return s
print("dG/dlogx numeric:", mp.nstr(dlogx, 25))
print("Li_{0,1}(x,y)   :", mp.nstr(Li01(x,y), 25))
print("dG/dlogy numeric:", mp.nstr(dlogy, 25))
print("Li_{1,0}(x,y)   :", mp.nstr(Li10(x,y), 25))

# Candidate DE C1: dG = G(b;z) dlog((b-a)/a) + G(a;z) dlog((z-b)/(b-a)), a=1/y, b=1/(xy)
a = 1/y; b = 1/(x*y); z = mp.mpf(1)
# dlog((b-a)/a): b-a = (1-x)/(xy), a=1/y -> ratio (1-x)/x ; dlog = -dlog(1-x)?? no: dlog((1-x)/x) = dlog(1-x)-dlog x
# coefficient of dlog x in dlog((1-x)/x): -1  => contrib G(b;z)*(-1)
# dlog((z-b)/(b-a)) = dlog(1-b) - dlog(b-a) = dlog((xy-1)/(xy)) - dlog((1-x)/(xy))
#   = dlog(1-xy) - dlog(1-x)   [constants cancel]
# coefficient of dlog x: 0 => contrib 0 wrt dlog x
Gb = mp.log(1-1/b)  # G(b;z)=log(1-z/b)... careful G(b;1)=log(1-1/b)
Ga = mp.log(1-1/a)
print("C1 predicts dG/dlogx = -G(b;z) =", mp.nstr(-Gb, 25))
print("actual dG/dlogx              =", mp.nstr(dlogx, 25))
