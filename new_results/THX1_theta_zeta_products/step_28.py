
import mpmath as mp
mp.mp.dps = 50
import pickle
from fractions import Fraction as Fr

# From the PROVED Dtilde formula, summing the geometric series:
#   D~(tau;m) = 2^{5/2} int_0^inf (gd x)^{2m} H(x) [(1+e^{-x/tau})^tau - 1] dx,  H(x) = e^{-x}(1+e^{-2x})^{-3/2}
def gd(x): return mp.atan(mp.sinh(x))
def H(x): return mp.e**(-x)*(1+mp.e**(-2*x))**mp.mpf('-1.5')

def Dtilde_GF2(tau, m):
    tau = mp.mpf(tau)
    def f(x):
        return gd(x)**(2*m) * H(x) * ((1+mp.e**(-x/tau))**tau - 1)
    return 2**mp.mpf('2.5')*mp.quad(f, [0, 2, 8, mp.inf])

Dtarget = pickle.load(open('Dtilde_relations.pkl','rb'))
def zval_target(n, m):
    coeffs = [Fr(c) for c in Dtarget[f"{n},{m}"]]
    if n==2: B=[mp.zeta(2)]
    elif n==3: B=[mp.zeta(3)]
    elif n==4: B=[mp.zeta(4)]
    elif n==5: B=[mp.zeta(5), mp.zeta(2)*mp.zeta(3)]
    elif n==6: B=[mp.zeta(6), mp.zeta(3)**2]
    elif n==7: B=[mp.zeta(7), mp.zeta(2)*mp.zeta(5), mp.zeta(2)**2*mp.zeta(3)]
    elif n==8: B=[mp.zeta(8), mp.zeta(3)*mp.zeta(5), mp.zeta(2)*mp.zeta(3)**2]
    elif n==9: B=[mp.zeta(9), mp.zeta(2)*mp.zeta(7), mp.zeta(2)**2*mp.zeta(5), mp.zeta(2)**3*mp.zeta(3), mp.zeta(3)**3]
    return sum(mp.mpf(c.numerator)/mp.mpf(c.denominator)*b for c,b in zip(coeffs,B))

for tau in [0.1, 0.05]:
    for m in [0, 1]:
        a = Dtilde_GF2(tau, m)
        tot = sum(zval_target(n,m)*mp.mpf(tau)**n for n in range(2,10) if n >= 2*m+2)
        print(f"tau={tau} m={m}: GF={mp.nstr(a,25)}  series={mp.nstr(tot,25)}  diff={mp.nstr(abs(a-tot),4)}")
