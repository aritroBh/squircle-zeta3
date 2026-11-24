import mpmath
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr
mp.dps = 70
def y_arc(x, p):
    # y>=0 solving |x|^p+y^p+|x+y|^p=1
    ax = abs(x)
    if ax**p >= 1: return mpf(0)
    lo, hi = mpf(0), mpf(1)
    for _ in range(160):
        mid = (lo+hi)/2
        if ax**p + mid**p + abs(x+mid)**p < 1: lo = mid
        else: hi = mid
    return (lo+hi)/2
def yp_arc(x, p):
    y = y_arc(x,p)
    sgn = lambda v: 1 if v>0 else -1
    num = sgn(x)*abs(x)**(p-1) + sgn(x+y)*abs(x+y)**(p-1)
    den = y**(p-1) + sgn(x+y)*abs(x+y)**(p-1)
    return -num/den
p = mpf(50)
xmax = mpf(2)**(-1/p)
Q1 = mpmath.quad(lambda x: sqrt(1+yp_arc(x,p)**2), [0, xmax/2, xmax])
Q2 = mpmath.quad(lambda x: sqrt(1+yp_arc(x,p)**2), [-xmax, -xmax/2, 0])
P = 2*(Q1+Q2)*1  # central symmetry: Q3=Q1, Q4=Q2? wait Q2 is half of 2nd quadrant...
print("Q1 =", nstr(Q1,30), " Q2 =", nstr(Q2,30))
print("2*(Q1+Q2) =", nstr(2*(Q1+Q2),30))
print("support P(50) = 7.912278464092482151719377673277517510033")
