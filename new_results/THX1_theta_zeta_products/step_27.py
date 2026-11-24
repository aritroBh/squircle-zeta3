
import mpmath as mp
mp.mp.dps = 50
import pickle
from fractions import Fraction as Fr

def gd(x):
    return mp.atan(mp.sinh(x))

def Dtilde_GF_y(tau, m):
    tau = mp.mpf(tau)
    def f(y):
        t = (1-y)/(1+y)
        prof = mp.sqrt((1+y)/2) * (1 + t**(1/(1-tau)))**((1-tau)/2) - 1
        return prof * mp.asin(y)**(2*m) / mp.sqrt(1-y**2)
    return 2*mp.sqrt(2)*mp.quad(f, [0, mp.mpf('0.5'), 1])

def Dtilde_GF_x(tau, m):
    tau = mp.mpf(tau); p = 1/(1-tau)
    def f(x):
        hp = (1+mp.e**(-2*x))**mp.mpf('-0.5') * (1+mp.e**(-2*p*x))**(1/(2*p))
        return (hp-1) * gd(x)**(2*m) * mp.sech(x)
    return 2*mp.sqrt(2)*mp.quad(f, [0, 2, 8, mp.inf])

for tau in [0.1, 0.05]:
    a = Dtilde_GF_y(tau, 0); b = Dtilde_GF_x(tau, 0)
    print(tau, "y-int:", mp.nstr(a,20), " x-int:", mp.nstr(b,20), " diff:", mp.nstr(abs(a-b),5))

Dtarget = pickle.load(open('Dtilde_relations.pkl','rb'))
def zval_target(n, m):
    coeffs = [Fr(c) for c in Dtarget[f"{n},{m}"]]
    if n==2: B=[mp.zeta(2)]
    if n==3: B=[mp.zeta(3)]
    if n==4: B=[mp.zeta(4)]
    if n==5: B=[mp.zeta(5), mp.zeta(2)*mp.zeta(3)]
    if n==6: B=[mp.zeta(6), mp.zeta(3)**2]
    if n==7: B=[mp.zeta(7), mp.zeta(2)*mp.zeta(5), mp.zeta(2)**2*mp.zeta(3)]
    return sum(mp.mpf(c.numerator)/mp.mpf(c.denominator)*b for c,b in zip(coeffs,B))
tot = sum(zval_target(n,0)*mp.mpf('0.1')**n for n in range(2,10))
print("expected D~(0.1;0) ~", mp.nstr(tot,20))
