import mpmath
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr, ellipe
mp.dps = 70
def y_arc(x, p):
    if x**p >= 1: return mpf(0)
    lo, hi = mpf(0), mpf(1)
    for _ in range(150):
        mid = (lo+hi)/2
        if x**p + mid**p + (x+mid)**p < 1: lo = mid
        else: hi = mid
    return (lo+hi)/2
MATS = [((1,0),(0,1)),((0,1),(1,0)),((-1,-1),(0,1)),((0,1),(-1,-1)),((-1,-1),(1,0)),((1,0),(-1,-1)),
        ((-1,0),(0,-1)),((0,-1),(-1,0)),((1,1),(0,-1)),((0,-1),(1,1)),((1,1),(-1,0)),((-1,0),(1,1))]
def h_arc_dir(c, s, p):
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
    return (f1+f2)/2
def h_hex(theta, p):
    c, s = cos(theta), sin(theta)
    best = -mpf(10)
    for (a,b),(d,e) in MATS:
        v = h_arc_dir(a*c+d*s, b*c+e*s, p)
        if v > best: best = v
    return best
Pfull = mpmath.quad(lambda th: h_hex(th, mpf(2)), [0, pi/2, pi, 3*pi/2, 2*pi])
print("P full =", nstr(Pfull,40), " ellipse:", nstr(4*ellipe(mpf(2)/3),40))
for th in [mpf('0.1'), mpf('0.3')]:
    print("th:", nstr(h_hex(th,mpf(2)),25), " th+pi/3:", nstr(h_hex(th+pi/3,mpf(2)),25))
