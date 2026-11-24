import mpmath, json
from mpmath import mp, mpf, sqrt, pi, log, exp, cos, sin, nstr, ellipe
mp.dps = 70
def y_arc(x, p):
    # solve x^p + y^p + (x+y)^p = 1 for y>=0 (x>=0)
    if x**p >= 1: return mpf(0)
    lo, hi = mpf(0), mpf(1)
    for _ in range(150):
        mid = (lo+hi)/2
        if x**p + mid**p + (x+mid)**p < 1: lo = mid
        else: hi = mid
    return (lo+hi)/2
def transforms(x, y):
    # 12 dihedral images of the arc point
    return [(x,y),(y,x),(-x-y,y),(y,-x-y),(-x-y,x),(x,-x-y),
            (-x,-y),(-y,-x),(x+y,-y),(-y,x+y),(x+y,-x),(-x,x+y)]
def h_hex(theta, p):
    c, s = cos(theta), sin(theta)
    xmax = mpf(2)**(-1/p)
    def hof(x):
        return x*c + y_arc(x,p)*s
    lo, hi = mpf(0), xmax
    gr = (sqrt(5)-1)/2
    x1, x2 = hi-gr*(hi-lo), lo+gr*(hi-lo)
    f1, f2 = hof(x1), hof(x2)
    for _ in range(120):
        if f1 < f2:
            lo, x1, f1 = x1, x2, f2
            x2 = lo+gr*(hi-lo); f2 = hof(x2)
        else:
            hi, x2, f2 = x2, x1, f1
            x1 = hi-gr*(hi-lo); f1 = hof(x1)
    xbest = (x1+x2)/2; ybest = y_arc(xbest,p)
    return max(x*c+y*s for x,y in transforms(xbest,ybest))
# validate at p=2: ellipse 2x^2+2xy+2y^2=1: P = 4 E(e^2=2/3)... semi axes 1/sqrt(3)? and 1
P2 = 12*mpmath.quad(lambda th: h_hex(th, mpf(2)), [0, pi/6])
print("P(p=2) =", nstr(P2,40))
# ellipse perimeter: a=1, b=1/sqrt(3): P = 4 a E(1-b^2/a^2) with E complete 2nd kind (parameter m)
e2 = 1-mpf(1)/3
print("ellipse  =", nstr(4*ellipe(e2),40))
# p -> inf limit: hexagon perimeter 3(sqrt2+1)
print("hex limit:", nstr(3*(sqrt(2)+1),30))
# h at theta=0 for p=2: should be max x = 2^{-1/2}
print("h(0,p=2):", nstr(h_hex(mpf(0), mpf(2)),20), " vs 2^{-1/2}:", nstr(mpf(2)**mpf('-0.5'),20))
