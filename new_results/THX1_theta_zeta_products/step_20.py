
import mpmath as mp
mp.mp.dps = 60

def Dtilde_theta(n, m):
    f = lambda ps: mp.cos(ps)*(-2*mp.log(mp.cos(ps)))**n * (mp.pi/2-2*ps)**(2*m)
    return 4*mp.sqrt(2)/mp.factorial(n) * mp.quad(f, [0, mp.pi/4])

def Dtilde_yint(n, m):
    f = lambda y: mp.log(2/(1+y))**n * mp.asin(y)**(2*m) / mp.sqrt(1-y)
    return 2/mp.factorial(n) * mp.quad(f, [0, 1])

import pickle
Dtarget = pickle.load(open('Dtilde_relations.pkl','rb'))
from fractions import Fraction as Fr
def zbasis(n):
    if n in (2,3,4): return [mp.zeta(n)]
    if n == 5: return [mp.zeta(5), mp.zeta(2)*mp.zeta(3)]
def Dtilde_target(n, m):
    coeffs = [Fr(c) for c in Dtarget[f"{n},{m}"]]
    return sum(mp.mpf(c.numerator)/mp.mpf(c.denominator)*b for c, b in zip(coeffs, zbasis(n)))

for (n,m) in [(2,0),(3,0),(4,0),(4,1),(5,0),(5,1)]:
    print((n,m), "theta:", mp.nstr(Dtilde_theta(n,m),20))
    print("     yint: ", mp.nstr(Dtilde_yint(n,m),20), " target:", mp.nstr(Dtilde_target(n,m),20))
