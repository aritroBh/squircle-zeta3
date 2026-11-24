import mpmath
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr
mp.dps = 70
def y_arc(x, p):
    if x**p >= 1: return mpf(0)
    lo, hi = mpf(0), mpf(1)
    for _ in range(160):
        mid = (lo+hi)/2
        if x**p + mid**p + (x+mid)**p < 1: lo = mid
        else: hi = mid
    return (lo+hi)/2
def yp_arc(x, p):
    y = y_arc(x,p)
    return -(x**(p-1)+(x+y)**(p-1))/(y**(p-1)+(x+y)**(p-1))
def P_arc(p):
    xmax = mpf(2)**(-1/p)
    return 4*mpmath.quad(lambda x: sqrt(1+yp_arc(x,p)**2), [0, xmax/2, xmax])
# The curve has 4-fold symmetry (swap x,y and central): perimeter = 4 * arc in first quadrant? 
# First quadrant arc from (2^{-1/p},0) to (0,2^{-1/p}) -- by swap symmetry the full curve =
# 2 * (arc + its swap) ... let me just compute full perimeter via support function already validated
print("P_arc(50) [1st quadrant arc *4]:", nstr(P_arc(mpf(50)),40))
# hexagon support function at p=inf: min_z3 max(|c1-z3|,|c2-z3|,|z3|)
def h_hex_inf(c1, c2):
    def F(z):
        return max(abs(c1-z), abs(c2-z), abs(z))
    lo, hi = -mpf(2), mpf(2)
    gr = (sqrt(5)-1)/2
    x1, x2 = hi-gr*(hi-lo), lo+gr*(hi-lo)
    f1, f2 = F(x1), F(x2)
    for _ in range(200):
        if f1 < f2:
            lo, x1, f1 = x1, x2, f2
            x2 = lo+gr*(hi-lo); f2 = F(x2)
        else:
            hi, x2, f2 = x2, x1, f1
            x1 = hi-gr*(hi-lo); f1 = F(x1)
    return (f1+f2)/2
Pinf = 4*mpmath.quad(lambda th: h_hex_inf(cos(th),sin(th)), [0, pi/4, pi/2])
print("P hexagon via support p=inf:", nstr(Pinf,40), " vs 2sqrt2+4:", nstr(2*sqrt(2)+4,30))
