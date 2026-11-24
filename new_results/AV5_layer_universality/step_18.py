import mpmath
from mpmath import mp, mpf, sqrt, zeta, log, exp, quad, tanh, nstr
mp.dps = 80
g1 = lambda s: log(1+exp(-2*s))
def Iab(a,b,splits):
    return quad(lambda s: s**a*g1(s)**b, splits)
for splits in [[0,1,3,8,20,45],[0,mpf('0.5'),1,2,4,8,16,32],[0,1,2,4,8,16,32,64]]:
    print("I(0,3):", nstr(Iab(0,3,splits),55))
    print("I(1,2):", nstr(Iab(1,2,splits),55))
# analytic pieces
I21 = mpf(7)/32*zeta(4)
Ig3 = -mpf(7)/16*zeta(4)
print("F3 from pieces:", nstr(Ig3 + Iab(0,3,[0,1,3,8,20,45])/6 + Iab(1,2,[0,1,3,8,20,45])/2 - I21/2, 55))
